% -------------------------------------------------------------------------
% PROJECT: Brain Tumor Edge Detection (Zotin et al. Reproduction)
% SCRIPT: main.m
% DATABASE: Task01_BrainTumour (MSD Structure)
% DESCRIPTION: Pipeline for T2 processing, Masking, BCET and GT loading.
% -------------------------------------------------------------------------

clear; clc; close all;

% --- STEP 1: Dynamic Path Configuration ---
% Determine the directory where this script resides (src folder)
script_dir = fileparts(mfilename('fullpath'));

% Navigate one level up to the project root, then into the data folder
data_dir = fullfile(script_dir, '..', 'data'); 
patient_id = 'BRATS_001.nii.gz';

image_path = fullfile(data_dir, 'imagesTr', patient_id);
label_path = fullfile(data_dir, 'labelsTr', patient_id);

% Verify file existence before loading to ensure pipeline robustness
if ~exist(image_path, 'file')
    error('Image file not found. Check if the dataset is extracted inside the correct "data" folder.');
end

% =========================================================================
% STEP 2: DATA LOADING AND SLICE EXTRACTION
% =========================================================================

disp('--- LOADING NIFTI VOLUMES ---');

% Read image metadata and volume (4D tensor: [X, Y, Z, Modality])
info_img = niftiinfo(image_path);
volume4D = niftiread(info_img);

% Read ground truth metadata and volume (3D tensor: [X, Y, Z])
info_lbl = niftiinfo(label_path);
volume_gt = niftiread(info_lbl);

% Extract the T2 modality channel (Index 4 for MSD Task01)
volume_t2 = double(volume4D(:, :, :, 4));

% Calculate the central axial slice index along the Z-axis
z_slices = size(volume_t2, 3);
mid_slice = round(z_slices / 2);

% Extract the 2D slices for both the image and the ground truth
slice_t2 = volume_t2(:, :, mid_slice);
slice_gt = double(volume_gt(:, :, mid_slice));

% Min-Max Normalization of the T2 slice to [0, 1] range
slice_norm = (slice_t2 - min(slice_t2(:))) / (max(slice_t2(:)) - min(slice_t2(:)));

disp(['Successfully extracted slice: ', num2str(mid_slice), ' out of ', num2str(z_slices)]);

% =========================================================================
% STEP 3: BRAIN MASKING AND BCET ENHANCEMENT
% =========================================================================

disp('--- APPLYING PRE-PROCESSING ---');

% 1. Brain Mask Generation (Otsu Method)
level = graythresh(slice_norm);
initial_mask = imbinarize(slice_norm, level);

% 2. Morphological cleaning of the mask
se = strel('disk', 5);
clean_mask = imopen(initial_mask, se);
clean_mask = imclose(clean_mask, se);

% 3. Application of the Balance Contrast Enhancement Technique (BCET)
enhanced_slice = apply_bcet(slice_norm, clean_mask);

disp('Pre-processing completed successfully.');

% --- Visualizzazione Intermedia ---
figure('Name', 'Pre-processing Results');
subplot(1, 3, 1);
imshow(slice_norm, []);
title('Original T2 Slice');

subplot(1, 3, 2);
imshow(clean_mask);
title('Morphological Brain Mask');

subplot(1, 3, 3);
imshow(enhanced_slice, []);
title('BCET Enhanced Output');

% --- Preliminary Visualization ---
% Min-Max Normalization of the T2 slice to [0, 1] range
slice_norm = (slice_t2 - min(slice_t2(:))) / (max(slice_t2(:)) - min(slice_t2(:)));

disp(['Successfully extracted slice: ', num2str(mid_slice), ' out of ', num2str(z_slices)]);

% --- Preliminary Visualization ---
% Call the custom function to display the extracted slices
display_preliminary_visualization(slice_norm, slice_gt);