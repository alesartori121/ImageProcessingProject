function [Pcd, Pnd, Pfa, FOM, Sens, Acc] = calculate_zotin_metrics_optimized(pred_edges, gt_mask)
% CALCULATE_ZOTIN_METRICS_OPTIMIZED Edge evaluation with Dilated Ground Truth.
%   Introduces a tolerance band (1 pixel radius) to account for human 
%   tracer variability and Canny's strict 1-pixel thickness.

    % 1. Extract the strict edges from the Ground Truth mask
    gt_edges_strict = logical(edge(gt_mask, 'canny'));
    pred_edges = logical(pred_edges);

    % --- INNOVATION: DILATED GROUND TRUTH (Tolerance Band) ---
    % Create a structural element (a disk of radius 1 pixel)
    se = strel('disk', 1); 
    % Dilate the strict edge to create a 3-pixel wide valid band
    gt_edges_dilated = imdilate(gt_edges_strict, se);

    % 2. Calculate Confusion Matrix (Using Dilated GT for True Positives)
    % A detection is "Correct" if it falls anywhere inside the tolerance band
    TP = sum(pred_edges(:) & gt_edges_dilated(:));  
    
    % False Positives are edges that are strictly OUTSIDE the tolerance band
    FP = sum(pred_edges(:) & ~gt_edges_dilated(:)); 
    
    % False Negatives are strict GT edges missed by our prediction
    FN = sum(~pred_edges(:) & gt_edges_strict(:)); 
    
    TN = sum(~pred_edges(:) & ~gt_edges_dilated(:));

    % 3. Calculate Probabilities and Basic Metrics
    Pcd = TP / (TP + FN);
    Sens = Pcd; 
    Pnd = FN / (TP + FN); 
    Pfa = FP / (FP + TN); 
    Acc = (TP + TN) / (TP + TN + FP + FN);

    % 4. Calculate Pratt's Figure of Merit (FOM) - Unchanged, still uses strict edges
    N_I = sum(gt_edges_strict(:));   
    N_A = sum(pred_edges(:)); 

    if N_A == 0 || N_I == 0
        FOM = 0.0; 
    else
        dist_map = bwdist(gt_edges_strict);
        d = dist_map(pred_edges);
        alpha = 1/9;
        FOM = (1 / max(N_I, N_A)) * sum(1 ./ (1 + alpha * (d.^2)));
    end
end