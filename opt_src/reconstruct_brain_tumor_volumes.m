function [brain_mask_3D, tumor_mask_3D, voxel_spacing, diagnostics] = reconstruct_brain_tumor_volumes(image_path, num_clusters, growth_iou_threshold)
    if nargin < 2 || isempty(num_clusters)
        num_clusters = 4;
    end
    if nargin < 3 || isempty(growth_iou_threshold)
        growth_iou_threshold = 0.15;
    end
    MIN_BRAIN_AREA = 50; 
    info = niftiinfo(image_path);
    volume4D = niftiread(info);
    volume_flair = double(volume4D(:, :, :, 1));
    voxel_spacing = info.PixelDimensions(1:3);
    [num_rows, num_cols, num_slices] = size(volume_flair);
    brain_voxels_all = volume_flair(volume_flair > 0);
    vol_min = min(brain_voxels_all);
    vol_max = max(brain_voxels_all);
    volume_norm = (volume_flair - vol_min) / (vol_max - vol_min);
    volume_norm(volume_norm < 0) = 0;
    volume_norm(volume_norm > 1) = 1;
    brain_mask_3D = false(num_rows, num_cols, num_slices);
    filtered_volume = zeros(num_rows, num_cols, num_slices);
    processed_slices = 0;
    for z = 1:num_slices
        try
            slice_norm = volume_norm(:, :, z);
            if max(slice_norm(:)) - min(slice_norm(:)) < eps
                continue;
            end
            filtered_slice = medfilt2(slice_norm, [3 3]);
            level = graythresh(filtered_slice);
            initial_mask = imbinarize(filtered_slice, level);
            clean_mask = bwareafilt(imfill(initial_mask, 'holes'), 1);
            if sum(clean_mask(:)) < MIN_BRAIN_AREA
                continue;
            end
            brain_mask_3D(:, :, z) = clean_mask;
            filtered_volume(:, :, z) = filtered_slice;
            processed_slices = processed_slices + 1;
        catch ME
            warning('Skipping slice %d (brain masking): %s', z, ME.message);
        end
    end
    enhanced_volume = apply_bcet(filtered_volume, brain_mask_3D);
    fcm_options = [2.0; 100; 1e-5; 0];
    fcm_input = enhanced_volume(brain_mask_3D);
    [~, U] = fcm(fcm_input, num_clusters, fcm_options);
    [~, max_U_idx] = max(U, [], 1);
    segmented_volume = zeros(size(enhanced_volume));
    segmented_volume(brain_mask_3D) = max_U_idx;
    hyper_threshold = quantile(brain_voxels_all, 0.95);
    hyperintense_area_per_slice = squeeze(sum(sum(volume_flair > hyper_threshold, 1), 2));
    [~, best_slice_idx] = max(hyperintense_area_per_slice);
    seed_slice_flair = volume_flair(:, :, best_slice_idx);
    seed_slice_norm = (seed_slice_flair - min(seed_slice_flair(:))) / ...
        (max(seed_slice_flair(:)) - min(seed_slice_flair(:)));
    seed_filtered = medfilt2(seed_slice_norm, [3 3]);
    seed_level = graythresh(seed_filtered);
    seed_initial_mask = imbinarize(seed_filtered, seed_level);
    seed_clean_mask = bwareafilt(imfill(seed_initial_mask, 'holes'), 1);
    seed_enhanced = apply_bcet(seed_filtered, seed_clean_mask);
    [~, seed_candidate_mask] = apply_fcm_clustering(seed_enhanced, seed_clean_mask, num_clusters);
    seed_mask_2D = isolate_tumor_mass(seed_candidate_mask);
    best_overlap = -1;
    tumor_cluster_idx = 1;
    for c = 1:num_clusters
        cluster_slice = (segmented_volume(:, :, best_slice_idx) == c);
        overlap = nnz(cluster_slice & seed_mask_2D);
        if overlap > best_overlap
            best_overlap = overlap;
            tumor_cluster_idx = c;
        end
    end
    AREA_STOP_FRACTION = 0.15;
    AREA_JUMP_FACTOR = 2.5;
    seed_area = nnz(seed_mask_2D);
    candidate_tumor_mask = false(num_rows, num_cols, num_slices);
    candidate_tumor_mask(:, :, best_slice_idx) = seed_mask_2D;
    prev_mask = seed_mask_2D;
    for z = (best_slice_idx - 1):-1:1
        [grown_mask, has_candidate] = grow_one_slice(segmented_volume(:, :, z), tumor_cluster_idx, prev_mask, growth_iou_threshold);
        grown_area = nnz(grown_mask);
        if ~has_candidate || grown_area < AREA_STOP_FRACTION * seed_area || grown_area > AREA_JUMP_FACTOR * nnz(prev_mask)
            break;
        end
        candidate_tumor_mask(:, :, z) = grown_mask;
        prev_mask = grown_mask;
    end
    prev_mask = seed_mask_2D;
    for z = (best_slice_idx + 1):num_slices
        [grown_mask, has_candidate] = grow_one_slice(segmented_volume(:, :, z), tumor_cluster_idx, prev_mask, growth_iou_threshold);
        grown_area = nnz(grown_mask);
        if ~has_candidate || grown_area < AREA_STOP_FRACTION * seed_area || grown_area > AREA_JUMP_FACTOR * nnz(prev_mask)
            break;
        end
        candidate_tumor_mask(:, :, z) = grown_mask;
        prev_mask = grown_mask;
    end
    tumor_mask_3D = imfill(candidate_tumor_mask, 'holes');
    diagnostics.num_slices = num_slices;
    diagnostics.processed_slices = processed_slices;
    diagnostics.best_slice_idx = best_slice_idx;
    diagnostics.raw_tumor_voxels = nnz(tumor_mask_3D);
    se = strel('sphere', 3);
    smoothed_tumor_mask = imopen(imclose(tumor_mask_3D, se), se);
    if nnz(smoothed_tumor_mask) > 0 || nnz(tumor_mask_3D) == 0
        tumor_mask_3D = smoothed_tumor_mask;
    end
    diagnostics.clean_tumor_voxels = nnz(tumor_mask_3D);
end

function [grown_mask, has_candidate] = grow_one_slice(cluster_volume_slice, tumor_cluster_idx, prev_mask, iou_threshold)
    cluster_slice = (cluster_volume_slice == tumor_cluster_idx);
    has_candidate = any(cluster_slice(:));
    grown_mask = false(size(cluster_slice));
    if ~has_candidate
        return;
    end
    prev_area = nnz(prev_mask);
    cc = bwconncomp(cluster_slice, 8);
    for ci = 1:cc.NumObjects
        comp_mask = false(size(cluster_slice));
        comp_mask(cc.PixelIdxList{ci}) = true;
        overlap_ratio = nnz(comp_mask & prev_mask) / prev_area;
        if overlap_ratio >= iou_threshold
            grown_mask = grown_mask | comp_mask;
        end
    end
    has_candidate = any(grown_mask(:));
end
