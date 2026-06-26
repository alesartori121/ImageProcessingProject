% -------------------------------------------------------------------------
% PROJECT: Brain Tumor Edge Detection
% SCRIPT: batch_processing.m
% DESCRIPTION: Executes the Zotin pipeline on ALL available BRATS patients,
% computes the metrics, and saves the results to a CSV file.
% -------------------------------------------------------------------------

clear; clc; close all;

disp('======================================================');
disp('   STARTING BATCH PROCESSING ON BRATS DATASET');
disp('======================================================');

% --- 1. SETUP PATHS ---
script_dir = fileparts(mfilename('fullpath'));
data_dir = fullfile(script_dir, '..', 'data'); 
images_dir = fullfile(data_dir, 'imagesTr');
labels_dir = fullfile(data_dir, 'labelsTr');

% Read all NIfTI files in the images directory
file_pattern = fullfile(images_dir, 'BRATS_*.nii.gz');
image_files = dir(file_pattern);
num_patients = length(image_files);

if num_patients == 0
    error('CRITICAL ERROR: No BRATS files found in the dataset folder.');
end

disp(['Found ', num2str(num_patients), ' patients to process.']);

% --- 2. PREALLOCATE DATA STORAGE ---
% Create an empty table to store the results
results_table = table('Size', [num_patients, 8], ...
    'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'PatientID', 'Pcd_Sens', 'Pnd', 'Pfa', 'Accuracy', 'FOM', 'ProcessingTime', 'NumClusters'});

% Initialize a Waitbar to track progress
h_wait = waitbar(0, 'Processing dataset, please wait...');

% --- 3. MAIN PROCESSING LOOP ---
for i = 1:num_patients
    tic; % Start timer for this patient
    
    patient_filename = image_files(i).name;
    results_table.PatientID(i) = patient_filename;
    
    img_path = fullfile(images_dir, patient_filename);
    lbl_path = fullfile(labels_dir, patient_filename);
    
    % Use try-catch to prevent a single bad file from stopping the whole loop
    try
        % -- Load Data --
        info_img = niftiinfo(img_path);
        vol_img = niftiread(info_img);
        info_lbl = niftiinfo(lbl_path);
        vol_lbl = double(niftiread(info_lbl));
        
        vol_t2 = double(vol_img(:, :, :, 4)); % T2 Modality
        mid_slice = round(size(vol_t2, 3) / 2);
        
        slice_t2 = vol_t2(:, :, mid_slice);
        slice_gt = vol_lbl(:, :, mid_slice);
        
        slice_norm = (slice_t2 - min(slice_t2(:))) / (max(slice_t2(:)) - min(slice_t2(:)));
        binary_gt = slice_gt > 0;
        
        % -- Pipeline --
        % NOTE ON FIDELITY: brain-masking (Otsu+morphology) below is an
        % addition vs. Zotin et al. (needed for BRATS' large background);
        % see src/main.m for the full rationale of every deviation.
        filtered_slice = medfilt2(slice_norm, [3 3]);
        level = graythresh(filtered_slice);
        initial_mask = imbinarize(filtered_slice, level);
        clean_mask = bwareafilt(imfill(initial_mask, 'holes'), 1);

        enhanced_slice = apply_bcet(filtered_slice, clean_mask);

        % c is chosen automatically from this range (Zotin et al. does not
        % specify c); see apply_fcm_clustering.m for the selection criteria.
        [~, candidate_tumor_mask, chosen_c] = apply_fcm_clustering(enhanced_slice, clean_mask, 3:5);

        final_tumor_mask = isolate_tumor_mass(candidate_tumor_mask);

        [tumor_edges, ~] = extract_tumor_edges(final_tumor_mask, slice_norm);

        % -- Evaluation --
        [Pcd, Pnd, Pfa, FOM, Sens, Acc] = calculate_zotin_metrics(tumor_edges, binary_gt);

        % -- Store Data --
        results_table.Pcd_Sens(i) = Sens; % Pcd is mathematically equal to Sensitivity
        results_table.Pnd(i) = Pnd;
        results_table.Pfa(i) = Pfa;
        results_table.Accuracy(i) = Acc;
        results_table.FOM(i) = FOM;
        results_table.NumClusters(i) = chosen_c;

    catch ME
        % If something fails (e.g., no tumor in the slice), log it and put NaNs
        warning('Error processing %s: %s', patient_filename, ME.message);
        results_table{i, 2:6} = NaN;
        results_table.NumClusters(i) = NaN;
    end
    
    results_table.ProcessingTime(i) = toc; % End timer
    
    % Update UI Waitbar
    waitbar(i / num_patients, h_wait, sprintf('Processing Patient %d of %d', i, num_patients));
end

close(h_wait);

% --- 4. SAVE AND SUMMARIZE RESULTS ---
disp('======================================================');
disp('   BATCH PROCESSING COMPLETE');
disp('======================================================');

% Save table to CSV
csv_filename = fullfile(script_dir, 'brats_statistical_results.csv');
writetable(results_table, csv_filename);
disp(['Raw data successfully saved to: ', csv_filename]);
disp(' ');

% Compute Mean and Standard Deviation (ignoring NaNs)
mean_sens = mean(results_table.Pcd_Sens, 'omitnan');
std_sens  = std(results_table.Pcd_Sens, 'omitnan');

mean_acc = mean(results_table.Accuracy, 'omitnan');
std_acc  = std(results_table.Accuracy, 'omitnan');

mean_fom = mean(results_table.FOM, 'omitnan');
std_fom  = std(results_table.FOM, 'omitnan');

% Print Summary
disp('--- OVERALL DATASET METRICS (MEAN ± STD) ---');
disp(['Sensitivity (Pcd) : ', num2str(mean_sens, '%.4f'), ' ± ', num2str(std_sens, '%.4f')]);
disp(['Accuracy          : ', num2str(mean_acc, '%.4f'),  ' ± ', num2str(std_acc, '%.4f')]);
disp(['Figure of Merit   : ', num2str(mean_fom, '%.4f'),  ' ± ', num2str(std_fom, '%.4f')]);
disp('--------------------------------------------');