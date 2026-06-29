
clear; clc; close all;

disp('======================================================');
disp('   STARTING BATCH PROCESSING ON BRATS DATASET');
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
    error('CRITICAL ERROR: No BRATS files found in the dataset folder.');
end
disp(['Found ', num2str(num_patients), ' patients to process.']);
% --- 2. PREALLOCATE DATA STORAGE ---
results_table = table('Size', [num_patients, 8], ...
    'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'PatientID', 'Pcd_Sens', 'Pnd', 'Pfa', 'Accuracy', 'FOM', 'ProcessingTime', 'NumClusters'});
% --- 3. MAIN PROCESSING LOOP ---
for i = 1:num_patients
    tic;
    patient_filename = image_files(i).name;
    results_table.PatientID(i) = patient_filename;
    img_path = fullfile(images_dir, patient_filename);
    lbl_path = fullfile(labels_dir, patient_filename);
    try
        % -- Load Data --
        info_img = niftiinfo(img_path);
        vol_img = niftiread(info_img);
        info_lbl = niftiinfo(lbl_path);
        vol_lbl = double(niftiread(info_lbl));
        vol_t2 = double(vol_img(:, :, :, 4));
        mid_slice = round(size(vol_t2, 3) / 2);
        slice_t2 = vol_t2(:, :, mid_slice);
        slice_gt = vol_lbl(:, :, mid_slice);
        slice_norm = (slice_t2 - min(slice_t2(:))) / (max(slice_t2(:)) - min(slice_t2(:)));
        binary_gt = slice_gt > 0;
        % -- Pipeline --
        filtered_slice = medfilt2(slice_norm, [3 3]);
        level = graythresh(filtered_slice);
        initial_mask = imbinarize(filtered_slice, level);
        clean_mask = bwareafilt(imfill(initial_mask, 'holes'), 1);
        enhanced_slice = apply_bcet(filtered_slice, clean_mask);
        num_clusters = 4;
        [~, candidate_tumor_mask] = apply_fcm_clustering(enhanced_slice, clean_mask, num_clusters);
        final_tumor_mask = isolate_tumor_mass(candidate_tumor_mask);
        [tumor_edges, ~] = extract_tumor_edges(final_tumor_mask, slice_norm);
        % -- Evaluation --
        [Pcd, Pnd, Pfa, FOM, Sens, Acc] = calculate_zotin_metrics(tumor_edges, binary_gt);
        % -- Store Data --
        results_table.Pcd_Sens(i) = Sens;
        results_table.Pnd(i) = Pnd;
        results_table.Pfa(i) = Pfa;
        results_table.Accuracy(i) = Acc;
        results_table.FOM(i) = FOM;
        results_table.NumClusters(i) = num_clusters;
    catch ME
        % If something fails (e.g., no tumor in the slice), log it and put NaNs
        warning('Error processing %s: %s', patient_filename, ME.message);
        results_table{i, 2:6} = NaN;
        results_table.NumClusters(i) = NaN;
    end
    results_table.ProcessingTime(i) = toc; % End timer
    if mod(i, 10) == 0 || i == num_patients
        fprintf('Processing Patient %d of %d\n', i, num_patients);
    end
end
disp('======================================================');
disp('   BATCH PROCESSING COMPLETE');
disp('======================================================');
csv_filename = fullfile(script_dir, 'brats_statistical_results.csv');
writetable(results_table, csv_filename);
disp(['Raw data successfully saved to: ', csv_filename]);
disp(' ');
mean_sens = mean(results_table.Pcd_Sens, 'omitnan');
std_sens  = std(results_table.Pcd_Sens, 'omitnan');
mean_acc = mean(results_table.Accuracy, 'omitnan');
std_acc  = std(results_table.Accuracy, 'omitnan');
mean_fom = mean(results_table.FOM, 'omitnan');
std_fom  = std(results_table.FOM, 'omitnan');
disp('--- OVERALL DATASET METRICS (MEAN ± STD) ---');
disp(['Sensitivity (Pcd) : ', num2str(mean_sens, '%.4f'), ' ± ', num2str(std_sens, '%.4f')]);
disp(['Accuracy          : ', num2str(mean_acc, '%.4f'),  ' ± ', num2str(std_acc, '%.4f')]);
disp(['Figure of Merit   : ', num2str(mean_fom, '%.4f'),  ' ± ', num2str(std_fom, '%.4f')]);
disp('--------------------------------------------');
