% -------------------------------------------------------------------------
% PROJECT: Brain Tumor Edge Detection (Optimized Pipeline)
% SCRIPT: main_optimized.m
% DESCRIPTION: Introduces FLAIR modality and Dynamic Slice Selection 
%              to overcome the topological limitations of the baseline.
% -------------------------------------------------------------------------
clear; clc; close all;

disp('======================================================');
disp('   RUNNING OPTIMIZED PIPELINE (FLAIR + DYNAMIC SLICE) ');
disp('======================================================');

% --- 1. DYNAMIC PATH CONFIGURATION ---
script_dir = fileparts(mfilename('fullpath'));
data_dir = fullfile(script_dir, '..', 'data');
patient_id = 'BRATS_001.nii.gz';

image_path = fullfile(data_dir, 'imagesTr', patient_id);
label_path = fullfile(data_dir, 'labelsTr', patient_id);
patient_name = erase(patient_id, '.nii.gz');
plots_dir = fullfile(script_dir, '..', 'plots', 'opt_src', patient_name);

if ~exist(image_path, 'file')
    error('CRITICAL ERROR: Image file not found in data directory.');
end

% --- 2. DATA LOADING & DYNAMIC SLICE EXTRACTION ---
disp('--- LOADING VOLUMES & FINDING BEST SLICE ---');
volume4D = niftiread(niftiinfo(image_path));
volume_gt = double(niftiread(niftiinfo(label_path)));

% INNOVATION 1: Dynamic Slice Selection (UNSUPERVISED)
% NOTE ON FIDELITY: the previous version of this script picked the slice
% with the most ground-truth tumor pixels -- that is label leakage, since
% the ground truth would not be available for a real, unlabeled scan. The
% score below only uses FLAIR intensities: it finds a single hyperintensity
% threshold across the WHOLE brain volume (95th percentile of brain voxels),
% then picks the axial slice containing the most voxels above it. Tumor/
% edema is hyperintense on FLAIR relative to normal tissue, so this is a
% plausible (if imperfect) proxy for "the slice most likely to contain the
% lesion" that could run on a real scan with no label.
% INNOVATION 2: Extract FLAIR Modality (Channel 1) instead of T2 (Channel 4)
volume_flair = double(volume4D(:, :, :, 1));

brain_voxels = volume_flair(volume_flair > 0); % BRATS volumes are pre-skull-stripped
hyper_threshold = quantile(brain_voxels, 0.95);
hyperintense_area_per_slice = squeeze(sum(sum(volume_flair > hyper_threshold, 1), 2));
[~, best_slice_idx] = max(hyperintense_area_per_slice);

disp(['Dynamic Selection: Maximum FLAIR hyperintense area found at slice Z = ', num2str(best_slice_idx)]);

% Extract the optimal slice
slice_flair = volume_flair(:, :, best_slice_idx);
slice_gt = volume_gt(:, :, best_slice_idx);

% Min-Max Normalization [0, 1] and Ground Truth Binarization
slice_norm = (slice_flair - min(slice_flair(:))) / (max(slice_flair(:)) - min(slice_flair(:)));
binary_gt = slice_gt > 0;

% --- 3. PRE-PROCESSING & ENHANCEMENT ---
% NOTE ON FIDELITY: brain-masking (Otsu+morphology) below is an addition vs.
% Zotin et al.; see src/main.m for the full rationale of every deviation.
disp('--- APPLYING ENHANCEMENT ---');
filtered_slice = medfilt2(slice_norm, [3 3]);

% Robust brain mask
level = graythresh(filtered_slice);
initial_mask = imbinarize(filtered_slice, level);
clean_mask = bwareafilt(imfill(initial_mask, 'holes'), 1);

% Balance Contrast Enhancement Technique (BCET)
enhanced_slice = apply_bcet(filtered_slice, clean_mask);

% --- 4. FCM CLUSTERING & TOPOLOGICAL FILTERING ---
% c is chosen automatically from this range (Zotin et al. does not specify
% c); see apply_fcm_clustering.m for the selection criteria.
disp('--- EXECUTING FCM CLUSTERING ---');
cluster_range = 3:5;
[segmented_slice, candidate_tumor_mask, chosen_c] = apply_fcm_clustering(enhanced_slice, clean_mask, cluster_range);
disp(['FCM: selected c = ', num2str(chosen_c), ' (best partition coefficient in range)']);

final_tumor_mask = isolate_tumor_mass(candidate_tumor_mask);

% --- 5. CANNY EDGE DETECTION ---
disp('--- EXTRACTING TUMOR EDGES ---');
[tumor_edges, clinical_overlay] = extract_tumor_edges(final_tumor_mask, slice_norm);

% --- 6. QUANTITATIVE METRICS ---
disp('--- CALCULATING METRICS ---');
[Pcd, Pnd, Pfa, FOM, Sens, Acc] = calculate_zotin_metrics_optimized(tumor_edges, binary_gt);

disp('--------------------------------------------------');
disp('   OPTIMIZED SLICE EVALUATION RESULTS (BRATS_001) ');
disp('--------------------------------------------------');
disp(['Pcd (Prob. of Correct Detection): ', num2str(Pcd, '%.4f')]);
disp(['Pnd (Prob. of Non-Detection):     ', num2str(Pnd, '%.4f')]);
disp(['Pfa (Prob. of False Alarm):       ', num2str(Pfa, '%.6f')]);
disp(['Sensitivity:                      ', num2str(Sens, '%.4f')]);
disp(['Accuracy:                         ', num2str(Acc, '%.4f')]);
disp(['FOM (Pratt''s Figure of Merit):    ', num2str(FOM, '%.4f')]);
disp('--------------------------------------------------');

% --- 7. VISUALIZATION RENDERING ---
disp('Rendering visualizations...');
display_preliminary_visualization(slice_norm, slice_gt, filtered_slice, ...
    clean_mask, enhanced_slice, chosen_c, segmented_slice, ...
    candidate_tumor_mask, final_tumor_mask, tumor_edges, clinical_overlay, binary_gt, plots_dir);
disp(['Figures saved to: ', plots_dir]);