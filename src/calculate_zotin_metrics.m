function [Pcd, Pnd, Pfa, FOM, Sens, Acc] = calculate_zotin_metrics(pred_edges, gt_mask)
% CALCULATE_ZOTIN_METRICS Computes Edge Detection metrics including Pratt's FOM.
%   NOTE ON FIDELITY: gt_mask is binary, so its boundary is traced with the
%   same morphological edge detector used in extract_tumor_edges.m (course
%   slides, leaf2.m), not with edge(...,'canny') -- consistent with how the
%   predicted edges are produced, so the comparison is apples-to-apples.
    se = strel('square', 3);
    gt_edges = imdilate(gt_mask, se) & ~gt_mask;
    pred_edges = logical(pred_edges);

    TP = sum(pred_edges(:) & gt_edges(:));
    FP = sum(pred_edges(:) & ~gt_edges(:));
    FN = sum(~pred_edges(:) & gt_edges(:));
    TN = sum(~pred_edges(:) & ~gt_edges(:));

    RECnt = TP + FN; % Reference edge pixel count (Zotin et al., denominator of Eq. 10-12)
    Pcd = TP / RECnt;
    Sens = Pcd;
    Pnd = FN / RECnt;
    Pfa = FP / RECnt; % Eq. 12: normalized by RECnt, NOT (FP+TN) -- this is not a standard FPR
    Acc = (TP + TN) / (TP + TN + FP + FN);

    N_I = sum(gt_edges(:));
    N_A = sum(pred_edges(:));

    if N_A == 0 || N_I == 0
        FOM = 0.0; 
    else
        dist_map = bwdist(gt_edges);
        d = dist_map(pred_edges);
        alpha = 1/9;
        FOM = (1 / max(N_I, N_A)) * sum(1 ./ (1 + alpha * (d.^2)));
    end
end