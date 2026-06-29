function [Pcd, Pnd, Pfa, FOM, Sens, Acc] = calculate_zotin_metrics_optimized(pred_edges, gt_mask)
    se_edge = strel('square', 3);
    gt_edges_strict = imdilate(gt_mask, se_edge) & ~gt_mask;
    pred_edges = logical(pred_edges);
    se = strel('disk', 1); 
    gt_edges_dilated = imdilate(gt_edges_strict, se);
    TP = sum(pred_edges(:) & gt_edges_dilated(:));
    FP = sum(pred_edges(:) & ~gt_edges_dilated(:));
    FN = sum(~pred_edges(:) & gt_edges_strict(:));
    TN = sum(~pred_edges(:) & ~gt_edges_dilated(:));
    RECnt = sum(gt_edges_strict(:));
    Pcd = TP / (TP + FN);
    Sens = Pcd;
    Pnd = FN / (TP + FN);
    Pfa = FP / RECnt;
    Acc = (TP + TN) / (TP + TN + FP + FN);
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