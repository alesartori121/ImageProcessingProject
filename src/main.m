% -------------------------------------------------------------------------
% PROJECT: Brain Tumor Edge Detection (Zotin et al. Reproduction)
% SCRIPT: main.m
% DESCRIPTION: Single-patient computation pipeline with visualization.
% -------------------------------------------------------------------------
clear; clc; close all;

% --- 1. DYNAMIC PATH CONFIGURATION ---
disp('--- CONFIGURING PATHS ---');
script_dir = fileparts(mfilename('fullpath'));
data_dir = fullfile(script_dir, '..', 'data');
patient_id = 'BRATS_001.nii.gz';
image_path = fullfile(data_dir, 'imagesTr', patient_id);
label_path = fullfile(data_dir, 'labelsTr', patient_id);
patient_name = erase(patient_id, '.nii.gz');
plots_dir = fullfile(script_dir, '..', 'plots', 'src', patient_name);

if ~exist(image_path, 'file')
    error('CRITICAL ERROR: Image file not found.');
end

% --- 2. DATA LOADING AND SLICE EXTRACTION ---
disp('--- LOADING NIFTI VOLUMES ---');
volume4D = niftiread(niftiinfo(image_path));
volume_gt = niftiread(niftiinfo(label_path));

% Extract T2 Modality (Channel 4) and central slice
volume_t2 = double(volume4D(:, :, :, 4));
mid_slice = round(size(volume_t2, 3) / 2);
slice_t2 = volume_t2(:, :, mid_slice);
slice_gt = double(volume_gt(:, :, mid_slice));

% Normalization [0, 1] and Ground Truth Binarization
slice_norm = (slice_t2 - min(slice_t2(:))) / (max(slice_t2(:)) - min(slice_t2(:)));
binary_gt = slice_gt > 0;

% --- 3. PRE-PROCESSING & ENHANCEMENT ---
% NOTE ON FIDELITY: Zotin et al. only specify grayscale conversion before
% median filtering + BCET (Fig. 1); there is no brain-masking prestep in the
% paper. The Otsu mask below is an addition needed for BRATS volumes, which
% have a large all-zero background outside the skull (Zotin's own images
% appear pre-cropped). It restricts BCET's statistics and FCM's clustering
% to brain tissue only.
disp('--- APPLYING ENHANCEMENT ---');
% Median filter to remove equipment noise (Zotin et al., Eq. 1)
filtered_slice = medfilt2(slice_norm, [3 3]);

% ADDITION vs paper: Otsu thresholding + morphology to mask out background
level = graythresh(filtered_slice);
initial_mask = imbinarize(filtered_slice, level);
clean_mask = bwareafilt(imfill(initial_mask, 'holes'), 1);

% Balance Contrast Enhancement Technique (Zotin et al., Eq. 2-6)
enhanced_slice = apply_bcet(filtered_slice, clean_mask);

% --- 4. FUZZY C-MEANS (FCM) CLUSTERING & ISOLATION ---
% NOTE ON FIDELITY: the paper (Eq. 7-9) does not specify the number of
% clusters c, nor how to pick the tumor cluster after FCM. apply_fcm_clustering
% now (a) chooses c automatically from cluster_range by partition coefficient,
% and (b) scores the brightest clusters by connected-component solidity
% instead of always taking the single brightest one, since on T2 the
% brightest cluster is often CSF/ventricles rather than the tumor.
disp('--- EXECUTING FCM CLUSTERING ---');
cluster_range = 3:5;
[segmented_slice, candidate_tumor_mask, chosen_c] = apply_fcm_clustering(enhanced_slice, clean_mask, cluster_range);
disp(['FCM: selected c = ', num2str(chosen_c), ' (best partition coefficient in range)']);

% Isolate the largest connected mass (removes smaller CSF artifacts)
final_tumor_mask = isolate_tumor_mass(candidate_tumor_mask);

% --- 5. CANNY EDGE DETECTION ---
% NOTE ON FIDELITY: extract_tumor_edges runs edge(...,'canny') on a BINARY
% mask, which degenerates to perimeter tracing (no gradient/hysteresis steps
% remain to apply -- see Gonzalez & Woods Ch.10 for the 4-step Canny
% algorithm, and Zotin et al. Sec 2.2, which applies Canny to the segmented
% multi-level image, not to a pre-isolated binary blob).
disp('--- EXTRACTING TUMOR EDGES ---');
[tumor_edges, clinical_overlay] = extract_tumor_edges(final_tumor_mask, slice_norm);

% --- 6. QUANTITATIVE METRICS ---
disp('--- CALCULATING METRICS ---');
[Pcd, Pnd, Pfa, FOM, Sens, Acc] = calculate_zotin_metrics(tumor_edges, binary_gt);
disp('--------------------------------------------------');
disp('   SINGLE SLICE EVALUATION RESULTS (BRATS_001)    ');
disp('--------------------------------------------------');
disp(['Pcd (Prob. of Correct Detection): ', num2str(Pcd, '%.4f')]);
disp(['Pnd (Prob. of Non-Detection):     ', num2str(Pnd, '%.4f')]);
disp(['Pfa (Prob. of False Alarm):       ', num2str(Pfa, '%.6f')]);
disp(['Sensitivity:                      ', num2str(Sens, '%.4f')]);
disp(['Accuracy:                         ', num2str(Acc, '%.4f')]);
disp(['FOM (Pratt''s Figure of Merit):    ', num2str(FOM, '%.4f')]);
disp('--------------------------------------------------');

% --- 7. VISUALIZATION RENDERING ---
disp('Computation complete. Rendering visualizations...');
display_preliminary_visualization(slice_norm, slice_gt, filtered_slice, ...
    clean_mask, enhanced_slice, chosen_c, segmented_slice, ...
    candidate_tumor_mask, final_tumor_mask, tumor_edges, clinical_overlay, binary_gt, plots_dir);
disp(['Figures saved to: ', plots_dir]);
