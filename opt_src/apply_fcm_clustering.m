function [segmented_slice, candidate_tumor_mask, chosen_c] = apply_fcm_clustering(enhanced_slice, clean_mask, num_clusters_range)
% APPLY_FCM_CLUSTERING Segments brain tissue with FCM and isolates the
% cluster most likely to be the tumor.
%   NOTE ON FIDELITY: Zotin et al. (Eq. 7-9) do not specify the number of
%   clusters c, nor how to pick the tumor cluster once FCM has run. This
%   function (a) picks c automatically from num_clusters_range by minimizing
%   the Xie-Beni cluster validity index (Xie & Beni, 1991), and (b) scores
%   the two brightest clusters by the solidity of their largest connected
%   component (a compact, massive blob beats thin/scattered CSF), instead of
%   always taking the single brightest cluster, which is often CSF/
%   ventricles rather than the tumor.
%   The simpler FCM partition coefficient (mean(U.^2)) was tried first, but
%   it is known to decrease monotonically with c (Bezdek, 1981): it always
%   selected the smallest c in the range regardless of the image, which is
%   not a real data-driven choice. Xie-Beni balances cluster compactness
%   against the separation between cluster centers and does not have that
%   bias, so it is used here instead.
    if nargin < 3 || isempty(num_clusters_range)
        num_clusters_range = 3:5;
    end
    fcm_options = [2.0; 100; 1e-5; 0];
    fcm_input = enhanced_slice(clean_mask);
    m = fcm_options(1);

    best_xb = Inf;
    centers = []; U = []; chosen_c = num_clusters_range(1);
    for c = num_clusters_range
        [c_centers, c_U] = fcm(fcm_input, c, fcm_options);
        dist2 = (c_centers - fcm_input').^2; % c x N squared distances (1-D intensity data)
        center_gaps = abs(c_centers - c_centers.');
        center_gaps(logical(eye(c))) = Inf;
        min_center_dist2 = min(center_gaps(:))^2;
        xb = sum((c_U.^m) .* dist2, 'all') / (numel(fcm_input) * min_center_dist2);
        if xb < best_xb
            best_xb = xb;
            centers = c_centers; U = c_U; chosen_c = c;
        end
    end

    [~, max_U_idx] = max(U, [], 1);
    segmented_slice = zeros(size(enhanced_slice));
    segmented_slice(clean_mask) = max_U_idx;

    % Tumor candidates: among the brightest clusters (tumor/edema is
    % hyperintense, but so is CSF -- this is the ambiguity to resolve).
    % CAUTION: one of the brightest clusters is very often just "most of the
    % normal brain tissue" (gray+white matter), which is itself large AND
    % solid -- so scoring by solidity*area alone can pick the whole brain
    % instead of the tumor. A focal lesion should not occupy most of the
    % brain cross-section, so any cluster whose largest blob exceeds
    % MAX_AREA_FRACTION of the brain mask is excluded from candidacy first,
    % regardless of how solid/bright it is.
    MAX_AREA_FRACTION = 0.40;
    brain_area = sum(clean_mask(:));
    [~, sort_idx] = sort(centers, 'descend');
    num_candidates = min(3, chosen_c);
    best_score = -Inf;
    candidate_tumor_mask = false(size(enhanced_slice));
    for k = 1:num_candidates
        cluster_mask = (segmented_slice == sort_idx(k));
        largest_blob = bwareafilt(cluster_mask, 1);
        area = sum(largest_blob(:));
        if area == 0 || area > MAX_AREA_FRACTION * brain_area
            continue; % empty, or too large to plausibly be a focal lesion
        end
        props = regionprops(largest_blob, 'Solidity');
        score = props(1).Solidity * area; % compact, massive blob > thin/scattered CSF
        if score > best_score
            best_score = score;
            candidate_tumor_mask = cluster_mask;
        end
    end
    if best_score == -Inf
        % Fallback: no cluster passed the plausibility check -- revert to
        % the single brightest cluster (Zotin et al.'s original heuristic).
        candidate_tumor_mask = (segmented_slice == sort_idx(1));
    end
end
