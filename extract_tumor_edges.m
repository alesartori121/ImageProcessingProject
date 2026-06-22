function tumor_edges = extract_tumor_edges(img_segmented, img_bcet, brain_mask, num_clusters)
% EXTRACT_TUMOR_EDGES - Finds the tumor cluster and applies Canny edge detection
%   tumor_edges = extract_tumor_edges(img_segmented, img_bcet, brain_mask, num_clusters)

% 1. Find the cluster with the highest mean intensity
mean_intensities = zeros(1, num_clusters);
for k = 1:num_clusters
    % Calculate mean intensity ONLY inside the brain mask
    region_pixels = img_bcet(img_segmented == k & brain_mask == 1);
    if isempty(region_pixels)
        mean_intensities(k) = 0;
    else
        mean_intensities(k) = mean(region_pixels);
    end
end

% The cluster with the maximum intensity is assumed to be the tumor/fluid
[~, tumor_cluster_id] = max(mean_intensities);

% 2. Create the raw binary mask for the tumor
tumor_mask = (img_segmented == tumor_cluster_id);

% 3. Morphological Cleaning
tumor_mask_clean = bwareaopen(tumor_mask, 100); % Remove small noise
tumor_mask_clean = bwareafilt(tumor_mask_clean, 1); % Keep largest connected component
tumor_mask_clean = imfill(tumor_mask_clean, 'holes'); % Fill internal holes

% 4. Apply Canny Edge Detector
tumor_edges = edge(tumor_mask_clean, 'canny');
end