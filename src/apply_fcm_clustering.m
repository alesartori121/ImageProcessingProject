function [segmented_slice, candidate_tumor_mask] = apply_fcm_clustering(enhanced_slice, clean_mask, num_clusters)
    fcm_options = [2.0; 100; 1e-5; 0]; 
    fcm_input = enhanced_slice(clean_mask);
    [centers, U] = fcm(fcm_input, num_clusters, fcm_options);
    [~, max_U_idx] = max(U, [], 1);
    
    segmented_slice = zeros(size(enhanced_slice));
    segmented_slice(clean_mask) = max_U_idx;
    
    [~, sort_idx] = sort(centers, 'ascend');
    brightest_cluster_label = sort_idx(end);
    candidate_tumor_mask = (segmented_slice == brightest_cluster_label);
end