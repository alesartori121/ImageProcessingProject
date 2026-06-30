function [segmented_slice, candidate_tumor_mask] = apply_fcm_clustering(enhanced_slice, clean_mask, num_clusters)
    if nargin < 3 || isempty(num_clusters)
        num_clusters = 4;
    end
    fcm_options = [2.0; 100; 1e-5; 0];
    fcm_input = enhanced_slice(clean_mask);

    [centers, U] = fcm(fcm_input, num_clusters, fcm_options);

    [~, max_U_idx] = max(U, [], 1);
    segmented_slice = zeros(size(enhanced_slice));
    segmented_slice(clean_mask) = max_U_idx;

    MAX_AREA_FRACTION = 0.40;
    brain_area = sum(clean_mask(:));
    [~, sort_idx] = sort(centers, 'descend');
    num_candidates = min(3, num_clusters);
    best_score = -Inf;
    candidate_tumor_mask = false(size(enhanced_slice));
    for k = 1:num_candidates
        cluster_mask = (segmented_slice == sort_idx(k));
        largest_blob = bwareafilt(cluster_mask, 1);
        area = sum(largest_blob(:));
        if area == 0 || area > MAX_AREA_FRACTION * brain_area
            continue;
        end
        props = regionprops(largest_blob, 'Solidity');
        score = props(1).Solidity * area;
        if score > best_score
            best_score = score;
            candidate_tumor_mask = cluster_mask;
        end
    end
    if best_score == -Inf
        candidate_tumor_mask = (segmented_slice == sort_idx(1));
    end
end
