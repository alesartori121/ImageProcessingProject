% -------------------------------------------------------------------------
% PROJECT: Brain Tumor Edge Detection (Zotin et al. Reproduction)
% -------------------------------------------------------------------------

clear; clc; close all;

% --- STEP 1: Data Loading and Slice Extraction ---
file_path = 'Task01_BrainTumour.tar'; 

info = niftiinfo(file_path);
volume4D = niftiread(info);

% Extract FLAIR modality (typically index 1)
volume_t2 = double(volume4D(:, :, :, 4));

% Extract the central axial slice
z_slices = size(volume_t2, 3);
mid_slice = round(z_slices / 2);
slice_2d = volume_t2(:, :, mid_slice);

% Normalize the slice to [0, 1] range for processing
slice_norm = (slice_2d - min(slice_2d(:))) / (max(slice_2d(:)) - min(slice_2d(:)));

% --- STEP 2: Brain Mask Generation (Otsu + Morphology) ---
% Compute global threshold and create initial binary mask
level = graythresh(slice_norm);
initial_mask = imbinarize(slice_norm, level);

% Apply morphological operations to clean the mask
se = strel('disk', 5);
clean_mask = imopen(initial_mask, se);
clean_mask = imclose(clean_mask, se);

% --- STEP 3: Contrast Enhancement (BCET) ---
% Call the custom BCET function passing the image and the mask
enhanced_slice = apply_bcet(slice_norm, clean_mask);

% --- VISUALIZATION ---
figure('Name', 'Pre-processing Pipeline Results', 'Position', [100, 100, 1200, 400]);

subplot(1, 3, 1);
imshow(slice_norm, []);
title('Original FLAIR Slice');

subplot(1, 3, 2);
imshow(clean_mask);
title('Brain Mask');

subplot(1, 3, 3);
imshow(enhanced_slice, []);
title('BCET Enhanced Slice');