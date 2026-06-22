function visualize_pipeline_results(img_orig, img_masked, img_filtered, img_bcet, img_segmented_disp, num_clusters, tumor_edges)
% VISUALIZE_PIPELINE_RESULTS - Plots all the steps of the edge detection pipeline

figure('Name', 'Phase 1 & 2: Full Pipeline Optimized', 'NumberTitle', 'off', 'Position', [100, 100, 1800, 800]);

subplot(2, 3, 1);
imshow(img_orig);
title('Original Image');

subplot(2, 3, 2);
imshow(img_masked);
title('Brain Masking');

subplot(2, 3, 3);
imshow(img_filtered);
title('Median Filter');

subplot(2, 3, 4);
imshow(img_bcet);
title('Enhanced (BCET)');

subplot(2, 3, 5);
imshow(img_segmented_disp);
colormap(gca, 'parula');
title(['FCM Segmented (k=', num2str(num_clusters), ')']);

subplot(2, 3, 6);
imshow(tumor_edges);
title('Edge Map (Canny)');
end