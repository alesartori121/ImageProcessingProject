function [segmented_img, cluster_centers] = fcm_segmentation(img_input, num_clusters)
% FCM_SEGMENTATION - Applies Fuzzy C-Means clustering to an image
%   [segmented_img, cluster_centers] = fcm_segmentation(img_input, num_clusters)
%
%   This function reshapes the 2D image into a 1D dataset, applies the FCM
%   algorithm to group pixels with similar intensities, and reconstructs 
%   the segmented image based on the maximum membership degree.

    % 1. Get original image dimensions for later reconstruction
    [rows, cols] = size(img_input);

    % 2. Reshape the 2D image matrix into a 1D column vector (N x 1)
    data_vector = img_input(:);

    % 3. Apply the Fuzzy C-Means algorithm
    options = [2.0; 100; 1e-5; 0]; % Default FCM options
    [cluster_centers, U] = fcm(data_vector, num_clusters, options);

    % 4. Hardening (Defuzzification)
    % Assign each pixel to the cluster it belongs to the *most*
    [~, max_idx] = max(U);

    % 5. Reconstruct the 2D segmented image
    segmented_img = reshape(max_idx, rows, cols);

end