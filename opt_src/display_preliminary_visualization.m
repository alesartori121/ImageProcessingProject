function display_preliminary_visualization(slice_norm, slice_gt, filtered_slice, clean_mask, enhanced_slice, num_clusters, segmented_slice, candidate_tumor_mask, final_tumor_mask, tumor_edges, clinical_overlay, binary_gt, output_dir, modality_label)
    if nargin < 13
        output_dir = '';
    end
    if nargin < 14 || isempty(modality_label)
        modality_label = 'T2';
    end
    save_figures = ~isempty(output_dir);
    if save_figures && ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % --- 1. Data Understanding ---
    fig1 = figure('Name', ['Data Understanding: ', modality_label, ' and Ground Truth'], 'Position', [50, 50, 1200, 450]);
    subplot(1, 2, 1);
    imshow(slice_norm, []);
    title(['Normalized ', modality_label, ' Slice']);

    subplot(1, 2, 2);
    imshow(slice_gt, []);
    title('Original Ground Truth Labels');
    if save_figures
        exportgraphics(fig1, fullfile(output_dir, '01_data_understanding.png'));
    end

    % --- 2. Enhancement Pipeline (Step 3) ---
    fig2 = figure('Name', 'Step 3: Pre-processing & Enhancement', 'Position', [100, 100, 1200, 450]);
    subplot(1, 3, 1);
    imshow(filtered_slice, []);
    title(['Median Filtered ', modality_label]);

    subplot(1, 3, 2);
    imshow(clean_mask);
    title('Robust Brain Mask');

    subplot(1, 3, 3);
    imshow(enhanced_slice, []);
    title('BCET Enhanced Output');
    if save_figures
        exportgraphics(fig2, fullfile(output_dir, '02_enhancement.png'));
    end

    % --- 3. Segmentation Pipeline (Step 4) ---
    fig3 = figure('Name', 'Step 4: FCM Segmentation & Isolation', 'Position', [150, 150, 1200, 450]);
    subplot(1, 3, 1);
    imshow(segmented_slice, []);
    colormap(jet(num_clusters + 1));
    colorbar;
    title(['FCM Classes (C = ', num2str(num_clusters), ')']);

    subplot(1, 3, 2);
    imshow(candidate_tumor_mask);
    title('Candidate Tumor Cluster (Tumor + CSF)');

    subplot(1, 3, 3);
    imshow(final_tumor_mask);
    title('Isolated Tumor Mass');
    if save_figures
        exportgraphics(fig3, fullfile(output_dir, '03_segmentation.png'));
    end

    % --- 4. Edge Detection Details (Step 5) ---
    fig4 = figure('Name', 'Step 5: Edge Detection Details', 'Position', [200, 200, 1200, 450]);
    subplot(1, 3, 1);
    imshow(final_tumor_mask);
    title('Isolated Tumor Mass');

    subplot(1, 3, 2);
    imshow(tumor_edges);
    title('Binary Edge Map');

    subplot(1, 3, 3);
    imshow(clinical_overlay);
    title('Clinical Overlay (Edges on T2)');
    if save_figures
        exportgraphics(fig4, fullfile(output_dir, '04_edge_detection.png'));
    end

    % --- 5. Final Clinical Validation (COMBINED OVERLAY) ---
    fig5 = figure('Name', 'Final Validation: Edge Detection vs Ground Truth', 'Position', [250, 250, 1200, 450]);

    % Subplot A: Our Algorithmic Prediction (Red)
    subplot(1, 3, 1);
    imshow(clinical_overlay);
    title('Our Prediction (Red)');

    se = strel('square', 3);
    gt_edges = imdilate(binary_gt, se) & ~binary_gt;
    gt_overlay = repmat(slice_norm, [1, 1, 3]);
    r_gt = gt_overlay(:,:,1); g_gt = gt_overlay(:,:,2); b_gt = gt_overlay(:,:,3);
    r_gt(gt_edges) = 0; g_gt(gt_edges) = 1; b_gt(gt_edges) = 0;
    gt_overlay(:,:,1)=r_gt; gt_overlay(:,:,2)=g_gt; gt_overlay(:,:,3)=b_gt;

    subplot(1, 3, 2);
    imshow(gt_overlay);
    title('Radiologist Ground Truth (Green)');

    combined_overlay = repmat(slice_norm, [1, 1, 3]);
    r_comb = combined_overlay(:,:,1);
    g_comb = combined_overlay(:,:,2);
    b_comb = combined_overlay(:,:,3);

    overlap_edges = tumor_edges & gt_edges;
    only_pred = tumor_edges & ~gt_edges;
    only_gt = gt_edges & ~tumor_edges;

    r_comb(only_pred) = 1; g_comb(only_pred) = 0; b_comb(only_pred) = 0;

    r_comb(only_gt) = 0; g_comb(only_gt) = 1; b_comb(only_gt) = 0;

    r_comb(overlap_edges) = 1; g_comb(overlap_edges) = 1; b_comb(overlap_edges) = 0;
    combined_overlay(:,:,1) = r_comb;
    combined_overlay(:,:,2) = g_comb;
    combined_overlay(:,:,3) = b_comb;

    subplot(1, 3, 3);
    imshow(combined_overlay);
    title('Combined Overlay (Yellow=Match)');
    if save_figures
        exportgraphics(fig5, fullfile(output_dir, '05_final_validation.png'));
    end
    drawnow;
end
