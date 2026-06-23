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

% --- 3. PRE-PROCESSING & ENHANCEMENT (STRICT ZOTIN PIPELINE) ---
disp('--- APPLYING ENHANCEMENT ---');
% Median filter to remove equipment noise
filtered_slice = medfilt2(slice_norm, [3 3]);

% Robust brain masking (Otsu thresholding + Morphology)
level = graythresh(filtered_slice);
initial_mask = imbinarize(filtered_slice, level);
clean_mask = bwareafilt(imfill(initial_mask, 'holes'), 1);

% Apply Balance Contrast Enhancement Technique
enhanced_slice = apply_bcet(filtered_slice, clean_mask);

% --- 4. FUZZY C-MEANS (FCM) CLUSTERING & ISOLATION ---
disp('--- EXECUTING FCM CLUSTERING ---');
num_clusters = 4;
[segmented_slice, candidate_tumor_mask] = apply_fcm_clustering(enhanced_slice, clean_mask, num_clusters);

% Isolate the largest connected mass (removes smaller CSF artifacts)
final_tumor_mask = isolate_tumor_mass(candidate_tumor_mask);

% --- 5. CANNY EDGE DETECTION ---
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
    clean_mask, enhanced_slice, num_clusters, segmented_slice, ...
    candidate_tumor_mask, final_tumor_mask, tumor_edges, clinical_overlay, binary_gt);
