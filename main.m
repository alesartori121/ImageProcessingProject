% =========================================================================
% Project: Edge detection in MRI brain tumor images
% Script: main.m
% Author: Sartori Alessandro
% Description: Phase 1 - Image Loading, Grayscale Conversion, and Precision
% =========================================================================

clc; clear; close all;

%% --- 0. PATH CONFIGURATION ---
% Define the absolute or relative paths of the dataset.
% Modify these strings with the actual path.
base_dir = 'BRATS2012'; % Main folder
img_dir = fullfile(base_dir, 'Images', 'BRATS_001'); % Example: first subfolder
lbl_dir = fullfile(base_dir, 'Labels', 'BRATS_001');

% Get a list of all PNG files in the image directory
% 'dir' returns an array of structures with file information
image_files = dir(fullfile(img_dir, '*.png'));

% Check if the directory is empty
if isempty(image_files)
    error('Zero images found in this folder: %s', img_dir);
end

% Automatically select the first file from the list
% image_files(1).name gets the name of the first file found
target_slice = image_files(1).name; 

img_path = fullfile(img_dir, target_slice);
lbl_path = fullfile(lbl_dir, target_slice); % Useful for Phase 3

%% --- 1. LOADING AND PRE-PROCESSING ---

% Check if the file exists before proceeding
if ~exist(img_path, 'file')
    error('File not found: %s. Check the path.', img_path);
end

img_orig = imread(img_path);

% Dimension check and grayscale conversion
if size(img_orig, 3) == 3
    img_gray = rgb2gray(img_orig);
    disp('RGB image detected: converted to grayscale.');
else
    img_gray = img_orig;
    disp('Image already in grayscale.');
end

% Conversion to double precision to preserve the fidelity of clinical data [0.0, 1.0]
img_double = im2double(img_gray);

%% --- 1.5. BRAIN MASKING (Pre-processing Enhancement) ---
% Isolate the brain by calling an external function
brain_mask = create_brain_mask(img_double);

% Apply the mask: keep brain tissue, force background strictly to 0
img_double_masked = img_double .* brain_mask;
disp('Brain masking applied successfully to remove background interference.');

%% --- 2. SPATIAL FILTERING (Image Enhancement) ---

% From the Paper (Zotin et al.): Application of the Median Filter
% Motivation: Reduction of "salt and pepper" noise while preserving edges
% compared to linear filters (e.g., mean or Gaussian).

kernel_size = [3 3]; % Typical size to balance smoothing and edge preservation
img_filtered = medfilt2(img_double, kernel_size);
disp('Median Filter [3x3] applied successfully.');

% --- Application of BCET (Balance Contrast Enhancement Technique) ---
% Enhances the contrast for better segmentation in the next phase
img_bcet = apply_bcet(img_filtered);
disp('BCET enhancement applied successfully.');

%% --- 3. EDGE MAP GENERATION (Phase 2) ---

% --- Step 1: Fuzzy C-Means (FCM) Segmentation ---
% We segment the image into 'k' clusters based on pixel intensities.
% Typically, for a brain MRI with a tumor, we look for:
% Background, Cerebrospinal Fluid, Brain Tissue, and Tumor.
num_clusters = 4;
[img_segmented, centers] = fcm_segmentation(img_bcet, num_clusters);
disp(['FCM Segmentation completed with ', num2str(num_clusters), ' clusters.']);

% The 'img_segmented' output contains values ​​from 1 to num_clusters.
% To display it correctly, we scale it to the colormap.
img_segmented_disp = mat2gray(img_segmented);

% --- Step 2: Tumor Extraction and Edge Detection (Canny) ---
% Call external function to extract tumor mask and apply Canny
tumor_edges = extract_tumor_edges(img_segmented, img_bcet, brain_mask, num_clusters);
disp('Canny edge detection applied successfully to the cleaned tumor mask.');

%% --- PRELIMINARY RESULTS VISUALIZATION ---
% Visualizzazione dei risultati
visualize_pipeline_results(img_orig, img_double_masked, img_filtered, img_bcet, img_segmented_disp, num_clusters, tumor_edges);