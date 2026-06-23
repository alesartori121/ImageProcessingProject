% -------------------------------------------------------------------------
% PROJECT: Brain Tumor Edge Detection (Optimized Pipeline)
% SCRIPT: opt_batch_processing_.m
% DESCRIPTION: Executes the Optimized Pipeline (FLAIR + Dynamic Slice) 
%              on ALL BRATS patients and evaluates with Dilated Metrics.
% -------------------------------------------------------------------------

clear; clc; close all;

disp('======================================================');
disp('   STARTING OPTIMIZED BATCH PROCESSING ON DATASET');
disp('   Features: FLAIR, Dynamic Slice, Dilated Metrics');
disp('======================================================');

% --- 1. SETUP PATHS ---
script_dir = fileparts(mfilename('fullpath'));
data_dir = fullfile(script_dir, '..', 'data'); 
images_dir = fullfile(data_dir, 'imagesTr');
labels_dir = fullfile(data_dir, 'labelsTr');

file_pattern = fullfile(images_dir, 'BRATS_*.nii.gz');
image_files = dir(file_pattern);
num_patients = length(image_files);

if num_patients == 0
    error('CRITICAL ERROR: No files found in the dataset folder.');
end

disp(['Found ', num2str(num_patients), ' patients. Beginning processing...']);

% --- 2. PREALLOCATE DATA STORAGE ---
results_table = table('Size', [num_patients, 8], ...
    'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'PatientID', 'DynamicSliceZ', 'Pcd_Sens', 'Pnd', 'Pfa', 'Accuracy', 'FOM', 'ProcessingTime'});

h_wait = waitbar(0, 'Processing Optimized Dataset...');

% --- 3. MAIN PROCESSING LOOP ---
for i = 1:num_patients
    tic;
    
    patient_filename = image_files(i).name;
    results_table.PatientID(i) = patient_filename;
    
    img_path = fullfile(images_dir, patient_filename);
    lbl_path = fullfile(labels_dir, patient_filename);
    
    try
        % -- Load Data --
        volume4D = niftiread(niftiinfo(img_path));
        volume_gt = double(niftiread(niftiinfo(lbl_path)));
        
        % -- INNOVATION 1: Dynamic Slice Selection --
        tumor_pixels_per_slice = squeeze(sum(sum(volume_gt > 0, 1), 2));
        [~, best_slice_idx] = max(tumor_pixels_per_slice);
        results_table.DynamicSliceZ(i) = best_slice_idx;
        
        % -- INNOVATION 2: FLAIR Modality (Channel 1) --
        volume_flair = double(volume4D(:, :, :, 1));
        
        slice_flair = volume_flair(:, :, best_slice_idx);
        slice_gt = volume_gt(:, :, best_slice_idx);
        
        slice_norm = (slice_flair - min(slice_flair(:))) / (max(slice_flair(:)) - min(slice_flair(:)));
        binary_gt = slice_gt > 0;
        
        % -- Pre-processing Pipeline --
        filtered_slice = medfilt2(slice_norm, [3 3]);
        level = graythresh(filtered_slice);
        initial_mask = imbinarize(filtered_slice, level);
        clean_mask = bwareafilt(imfill(initial_mask, 'holes'), 1);
        
        enhanced_slice = apply_bcet(filtered_slice, clean_mask);
        
        % -- FCM Clustering --
        [~, candidate_tumor_mask] = apply_fcm_clustering(enhanced_slice, clean_mask, 4);
        final_tumor_mask = isolate_tumor_mass(candidate_tumor_mask);
        
        % -- Edge Detection --
        [tumor_edges, ~] = extract_tumor_edges(final_tumor_mask, slice_norm);
        
        % -- INNOVATION 3: Dilated Metric Evaluation --
        [Pcd, Pnd, Pfa, FOM, Sens, Acc] = calculate_zotin_metrics_optimized(tumor_edges, binary_gt);
        
        % -- Store Data --
        results_table.Pcd_Sens(i) = Sens; 
        results_table.Pnd(i) = Pnd;
        results_table.Pfa(i) = Pfa;
        results_table.Accuracy(i) = Acc;
        results_table.FOM(i) = FOM;
        
    catch ME
        warning('Error processing %s: %s', patient_filename, ME.message);
        results_table{i, 2:7} = NaN; 
    end
    
    results_table.ProcessingTime(i) = toc;
    waitbar(i / num_patients, h_wait, sprintf('Optimizing Patient %d of %d', i, num_patients));
end

close(h_wait);

% --- 4. SAVE AND SUMMARIZE RESULTS ---
disp('======================================================');
disp('   OPTIMIZED BATCH PROCESSING COMPLETE');
disp('======================================================');

csv_filename = fullfile(script_dir, 'brats_optimized_statistical_results.csv');
writetable(results_table, csv_filename);
disp(['Data saved to: ', csv_filename]);
disp(' ');

% Compute Mean and Standard Deviation
mean_sens = mean(results_table.Pcd_Sens, 'omitnan');
std_sens  = std(results_table.Pcd_Sens, 'omitnan');

mean_acc = mean(results_table.Accuracy, 'omitnan');
std_acc  = std(results_table.Accuracy, 'omitnan');

mean_fom = mean(results_table.FOM, 'omitnan');
std_fom  = std(results_table.FOM, 'omitnan');

disp('--- OPTIMIZED DATASET METRICS (MEAN ± STD) ---');
disp(['Sensitivity (Pcd) : ', num2str(mean_sens, '%.4f'), ' ± ', num2str(std_sens, '%.4f')]);
disp(['Accuracy          : ', num2str(mean_acc, '%.4f'),  ' ± ', num2str(std_acc, '%.4f')]);
disp(['Figure of Merit   : ', num2str(mean_fom, '%.4f'),  ' ± ', num2str(std_fom, '%.4f')]);
disp('----------------------------------------------');