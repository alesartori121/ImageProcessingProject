function [segmented_slice, candidate_tumor_mask] = apply_fcm_clustering(enhanced_slice, clean_mask, num_clusters)
% APPLY_FCM_CLUSTERING Segments brain tissue with FCM and isolates the
% cluster most likely to be the tumor.
%   NOTE ON FIDELITY: Zotin et al. (Eq. 7-9) do not specify the number of
%   clusters c, nor how to pick the tumor cluster once FCM has run. Gonzalez
%   & Woods (Ch. 10.5, clustering-based segmentation) treats the number of
%   clusters as a value that must be specified directly -- "the important
%   issue is the value selected for k ... multiple passes are rarely used"
%   -- rather than estimated automatically, so a fixed c=4 is used here.
%   The tumor cluster is then chosen among up to the three brightest
%   clusters by the solidity of their largest connected component (a
%   compact, massive blob beats thin/scattered CSF), instead of always
%   taking the single brightest cluster, which is often CSF/ventricles
%   rather than the tumor.
    if nargin < 3 || isempty(num_clusters)
        num_clusters = 4;
    end
    fcm_options = [2.0; 100; 1e-5; 0];
    fcm_input = enhanced_slice(clean_mask);

    [centers, U] = fcm(fcm_input, num_clusters, fcm_options);

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
    num_candidates = min(3, num_clusters);
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
