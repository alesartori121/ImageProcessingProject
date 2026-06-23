function [Pcd, Pnd, Pfa, FOM, Sens, Acc] = calculate_zotin_metrics(pred_edges, gt_mask)
% CALCULATE_ZOTIN_METRICS Computes Edge Detection metrics including Pratt's FOM.
    gt_edges = logical(edge(gt_mask, 'canny'));
    pred_edges = logical(pred_edges);

    TP = sum(pred_edges(:) & gt_edges(:));
    FP = sum(pred_edges(:) & ~gt_edges(:));
    FN = sum(~pred_edges(:) & gt_edges(:));
    TN = sum(~pred_edges(:) & ~gt_edges(:));

    Pcd = TP / (TP + FN);
    Sens = Pcd; 
    Pnd = FN / (TP + FN); 
    Pfa = FP / (FP + TN); 
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