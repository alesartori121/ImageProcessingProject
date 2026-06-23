function display_preliminary_visualization(slice_norm, slice_gt, filtered_slice, clean_mask, enhanced_slice, num_clusters, segmented_slice, candidate_tumor_mask, final_tumor_mask, tumor_edges, clinical_overlay, binary_gt)
% DISPLAY_PRELIMINARY_VISUALIZATION Renders all pipeline visualizations.
%   This function acts as the graphical front-end of the application, 
%   generating all the intermediate and final figures required for the 
%   scientific report, keeping the main script clean from UI logic.

    % --- 1. Data Understanding ---
    figure('Name', 'Data Understanding: T2 and Ground Truth', 'Position', [50, 50, 1000, 400]);
    subplot(1, 2, 1);
    imshow(slice_norm, []);
    title('Normalized T2 Slice');

    subplot(1, 2, 2);
    imshow(slice_gt, []);
    title('Original Ground Truth Labels');

    % --- 2. Enhancement Pipeline (Step 3) ---
    figure('Name', 'Step 3: Pre-processing & Enhancement', 'Position', [100, 100, 1200, 400]);
    subplot(1, 3, 1);
    imshow(filtered_slice, []);
    title('Median Filtered T2');

    subplot(1, 3, 2);
    imshow(clean_mask);
    title('Robust Brain Mask');

    subplot(1, 3, 3);
    imshow(enhanced_slice, []);
    title('BCET Enhanced Output');

    % --- 3. Segmentation Pipeline (Step 4) ---
    figure('Name', 'Step 4: FCM Segmentation & Isolation', 'Position', [150, 150, 1200, 400]);
    subplot(1, 3, 1);
    imshow(segmented_slice, []);
    colormap(jet(num_clusters + 1));
    colorbar;
    title(['FCM Classes (C = ', num2str(num_clusters), ')']);

    subplot(1, 3, 2);
    imshow(candidate_tumor_mask);
    title('Brightest Cluster (Tumor + CSF)');

    subplot(1, 3, 3);
    imshow(final_tumor_mask);
    title('Isolated Tumor Mass');

    % --- 4. Edge Detection Details (Step 5) ---
    figure('Name', 'Step 5: Edge Detection Details', 'Position', [200, 200, 1200, 450]);
    subplot(1, 3, 1);
    imshow(final_tumor_mask);
    title('Isolated Tumor Mass');

    subplot(1, 3, 2);
    imshow(tumor_edges);
    title('Binary Canny Edge Map');

    subplot(1, 3, 3);
    imshow(clinical_overlay);
    title('Clinical Overlay (Edges on T2)');

    % --- 5. Final Clinical Validation (COMBINED OVERLAY) ---
    figure('Name', 'Final Validation: Edge Detection vs Ground Truth', 'Position', [250, 250, 1400, 450]);
    
    % Subplot A: Our Algorithmic Prediction (Red)
    subplot(1, 3, 1);
    imshow(clinical_overlay);
    title('Our Prediction (Red)');

    % Subplot B: Ground Truth (Green)
    gt_edges = edge(binary_gt, 'canny');
    gt_overlay = repmat(slice_norm, [1, 1, 3]);
    r_gt = gt_overlay(:,:,1); g_gt = gt_overlay(:,:,2); b_gt = gt_overlay(:,:,3);
    r_gt(gt_edges) = 0; g_gt(gt_edges) = 1; b_gt(gt_edges) = 0;
    gt_overlay(:,:,1)=r_gt; gt_overlay(:,:,2)=g_gt; gt_overlay(:,:,3)=b_gt;
    
    subplot(1, 3, 2);
    imshow(gt_overlay);
    title('Radiologist Ground Truth (Green)');

    % Subplot C: Combined Mathematical Overlay (Yellow/Red/Green)
    combined_overlay = repmat(slice_norm, [1, 1, 3]);
    r_comb = combined_overlay(:,:,1); 
    g_comb = combined_overlay(:,:,2); 
    b_comb = combined_overlay(:,:,3);
    
    % Logic arrays for intersections and differences
    overlap_edges = tumor_edges & gt_edges;
    only_pred = tumor_edges & ~gt_edges;
    only_gt = gt_edges & ~tumor_edges;
    
    % Apply Red where ONLY the algorithm found an edge
    r_comb(only_pred) = 1; g_comb(only_pred) = 0; b_comb(only_pred) = 0;
    
    % Apply Green where ONLY the radiologist found an edge
    r_comb(only_gt) = 0; g_comb(only_gt) = 1; b_comb(only_gt) = 0;
    
    % Apply Yellow (Red+Green) where BOTH agree perfectly
    r_comb(overlap_edges) = 1; g_comb(overlap_edges) = 1; b_comb(overlap_edges) = 0;
    
    % Reconstruct the final image
    combined_overlay(:,:,1) = r_comb; 
    combined_overlay(:,:,2) = g_comb; 
    combined_overlay(:,:,3) = b_comb;

    subplot(1, 3, 3);
    imshow(combined_overlay);
    title('Combined Overlay (Yellow=Match)');
    
    % Force MATLAB to draw all figures immediately
    drawnow;
end