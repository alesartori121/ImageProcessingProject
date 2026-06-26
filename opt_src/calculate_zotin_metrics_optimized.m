function [Pcd, Pnd, Pfa, FOM, Sens, Acc] = calculate_zotin_metrics_optimized(pred_edges, gt_mask)
% CALCULATE_ZOTIN_METRICS_OPTIMIZED Edge evaluation with Dilated Ground Truth.
%   Introduces a tolerance band (1 pixel radius) to account for human
%   tracer variability and the strict 1-pixel thickness of the boundary.
%   NOTE ON FIDELITY: gt_mask is binary, so its boundary is traced with the
%   same morphological edge detector used in extract_tumor_edges.m (course
%   slides, leaf2.m), not with edge(...,'canny').

    % 1. Extract the strict edges from the Ground Truth mask
    se_edge = strel('square', 3);
    gt_edges_strict = imdilate(gt_mask, se_edge) & ~gt_mask;
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
    RECnt = sum(gt_edges_strict(:)); % Eq. 12 denominator (Zotin et al.): strict reference edge count
    Pcd = TP / (TP + FN); % Recall-style ratio; kept as TP/(TP+FN) rather than TP/RECnt because
                          % TP here is counted over predicted pixels matched within the dilated
                          % tolerance band, so it is not bounded by RECnt the way the strict
                          % metric is (see calculate_zotin_metrics.m for the strict version).
    Sens = Pcd;
    Pnd = FN / (TP + FN);
    Pfa = FP / RECnt; % Eq. 12: normalized by RECnt, NOT (FP+TN) -- this is not a standard FPR
    Acc = (TP + TN) / (TP + TN + FP + FN);

    % 4. Calculate Pratt's Figure of Merit (FOM) - Unchanged, still uses strict edges
    N_I = RECnt;
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