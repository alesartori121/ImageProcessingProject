% -------------------------------------------------------------------------
% PROJECT: Brain Tumor Edge Detection (Optimized Pipeline)
% SCRIPT: opt_batch_3d_reconstruction.m
% DESCRIPTION: Runs the 3-D brain+tumor reconstruction
%              (reconstruct_brain_tumor_volumes + render_brain_tumor_3d)
%              over every patient in the dataset, in parallel via parfor.
%              For each patient, saves a PNG preview and an interactive
%              MATLAB .fig (rotatable when reopened) under
%              plots/opt_src/<patient>/06_3d_reconstruction.{png,fig}.
% -------------------------------------------------------------------------
clear; clc; close all;

disp('======================================================');
disp('   BATCH 3-D BRAIN+TUMOR RECONSTRUCTION (ALL PATIENTS)');
disp('======================================================');

script_dir = fileparts(mfilename('fullpath'));
data_dir = fullfile(script_dir, '..', 'data');
images_dir = fullfile(data_dir, 'imagesTr');
plots_root = fullfile(script_dir, '..', 'plots', 'opt_src');
if ~exist(plots_root, 'dir')
    mkdir(plots_root);
end

file_pattern = fullfile(images_dir, 'BRATS_*.nii.gz');
image_files = dir(file_pattern);
num_patients = length(image_files);
if num_patients == 0
    error('CRITICAL ERROR: No files found in the dataset folder.');
end
disp(['Found ', num2str(num_patients), ' patients. Beginning parallel 3-D reconstruction...']);

num_clusters = 4;
results = cell(num_patients, 1);

parfor i = 1:num_patients
    patient_filename = image_files(i).name;
    patient_name = erase(patient_filename, '.nii.gz');
    img_path = fullfile(images_dir, patient_filename);
    patient_plots_dir = fullfile(plots_root, patient_name);
    if ~exist(patient_plots_dir, 'dir')
        mkdir(patient_plots_dir);
    end

    row = struct('PatientID', string(patient_filename), 'TotalSlices', NaN, ...
                  'ProcessedSlices', NaN, 'RawTumorVoxels', NaN, 'CleanTumorVoxels', NaN, ...
                  'ProcessingTime', NaN, 'ErrorMessage', "");
    t0 = tic;
    try
        [brain_mask_3D, tumor_mask_3D, voxel_spacing, diagnostics] = ...
            reconstruct_brain_tumor_volumes(img_path, num_clusters);

        fig = render_brain_tumor_3d(brain_mask_3D, tumor_mask_3D, voxel_spacing, patient_name, 'off');
        exportgraphics(fig, fullfile(patient_plots_dir, '06_3d_reconstruction.png'));
        savefig(fig, fullfile(patient_plots_dir, '06_3d_reconstruction.fig'));
        close(fig);

        row.TotalSlices = diagnostics.num_slices;
        row.ProcessedSlices = diagnostics.processed_slices;
        row.RawTumorVoxels = diagnostics.raw_tumor_voxels;
        row.CleanTumorVoxels = diagnostics.clean_tumor_voxels;
    catch ME
        row.ErrorMessage = string(ME.message);
        warning('Error reconstructing %s: %s', patient_filename, ME.message);
    end
    row.ProcessingTime = toc(t0);
    results{i} = row;
    fprintf('[%d/%d] %s done in %.1fs\n', i, num_patients, patient_filename, row.ProcessingTime);
end

% --- SAVE SUMMARY ---
disp('======================================================');
disp('   BATCH 3-D RECONSTRUCTION COMPLETE');
disp('======================================================');
results_table = struct2table([results{:}]);
csv_path = fullfile(script_dir, 'brats_3d_reconstruction_summary.csv');
writetable(results_table, csv_path);
disp(['Summary saved to: ', csv_path]);

num_failed = sum(results_table.ErrorMessage ~= "");
disp([num2str(num_patients - num_failed), ' / ', num2str(num_patients), ' patients reconstructed successfully.']);
if num_failed > 0
    disp([num2str(num_failed), ' patients failed (see ErrorMessage column in the summary CSV).']);
end
