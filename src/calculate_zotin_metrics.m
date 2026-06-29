function [Pcd, Pnd, Pfa, FOM, Sens, Acc] = calculate_zotin_metrics(pred_edges, gt_mask)
    se = strel('square', 3);
    gt_edges = imdilate(gt_mask, se) & ~gt_mask;
    pred_edges = logical(pred_edges);
    TP = sum(pred_edges(:) & gt_edges(:));
    FP = sum(pred_edges(:) & ~gt_edges(:));
    FN = sum(~pred_edges(:) & gt_edges(:));
    TN = sum(~pred_edges(:) & ~gt_edges(:));
    RECnt = TP + FN;
    Pcd = TP / RECnt;
    Sens = Pcd;
    Pnd = FN / RECnt;
    Pfa = FP / RECnt;
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