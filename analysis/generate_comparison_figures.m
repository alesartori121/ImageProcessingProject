% -------------------------------------------------------------------------
% PROJECT: Brain Tumor Edge Detection -- Ablation Study
% SCRIPT: generate_comparison_figures.m
% DESCRIPTION: Regenerates and saves the pipeline visualization figures for
% representative patients (best/median/worst FOM) of each pipeline, for use
% in the report -- qualitative figures no longer depend on whichever single
% patient happens to be hardcoded in main.m/opt_main.m.
% -------------------------------------------------------------------------
clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
project_dir = fullfile(script_dir, '..');
data_dir = fullfile(project_dir, 'data');

src_csv = fullfile(project_dir, 'src', 'brats_statistical_results.csv');
opt_csv = fullfile(project_dir, 'opt_src', 'brats_optimized_statistical_results.csv');

if ~exist(src_csv, 'file') || ~exist(opt_csv, 'file')
    error('CRITICAL ERROR: run src/batch_processing.m and opt_src/opt_batch_processing.m first.');
end

%% --- Baseline pipeline (src) ---
results_src = readtable(src_csv);
picks_src = pick_representative_cases(results_src);
src_dir = fullfile(project_dir, 'src');
addpath(src_dir);
for i = 1:numel(picks_src)
    pick = picks_src(i);
    fprintf('[src] %s case: %s (FOM=%.4f)\n', pick.label, pick.PatientID, pick.FOM);
    out_dir = fullfile(project_dir, 'plots', 'src', [pick.label, '_', erase(pick.PatientID, '.nii.gz')]);
    process_baseline_patient(data_dir, pick.PatientID, out_dir);
    close all;
end
rmpath(src_dir);

%% --- Optimized pipeline (opt_src) ---
results_opt = readtable(opt_csv);
picks_opt = pick_representative_cases(results_opt);
opt_dir = fullfile(project_dir, 'opt_src');
addpath(opt_dir);
for i = 1:numel(picks_opt)
    pick = picks_opt(i);
    fprintf('[opt_src] %s case: %s (FOM=%.4f)\n', pick.label, pick.PatientID, pick.FOM);
    out_dir = fullfile(project_dir, 'plots', 'opt_src', [pick.label, '_', erase(pick.PatientID, '.nii.gz')]);
    process_optimized_patient(data_dir, pick.PatientID, out_dir);
    close all;
end
rmpath(opt_dir);

disp('All comparison figures generated under plots/src and plots/opt_src.');

%% --- Local functions ---

function picks = pick_representative_cases(results_table)
% Picks the best, median, and worst patient by FOM (NaN rows excluded).
    valid = ~isnan(results_table.FOM);
    t = results_table(valid, :);
    [~, idx_best] = max(t.FOM);
    [~, idx_worst] = min(t.FOM);
    [~, idx_median] = min(abs(t.FOM - median(t.FOM)));

    picks(1) = struct('label', 'best',   'PatientID', t.PatientID(idx_best),   'FOM', t.FOM(idx_best));
    picks(2) = struct('label', 'median', 'PatientID', t.PatientID(idx_median), 'FOM', t.FOM(idx_median));
    picks(3) = struct('label', 'worst',  'PatientID', t.PatientID(idx_worst),  'FOM', t.FOM(idx_worst));
end

function process_baseline_patient(data_dir, patient_id, plots_dir)
% Mirrors src/main.m's pipeline for a single patient, saving figures only.
    image_path = fullfile(data_dir, 'imagesTr', patient_id);
    label_path = fullfile(data_dir, 'labelsTr', patient_id);

    volume4D = niftiread(niftiinfo(image_path));
    volume_gt = niftiread(niftiinfo(label_path));
    volume_t2 = double(volume4D(:, :, :, 4));
    mid_slice = round(size(volume_t2, 3) / 2);
    slice_t2 = volume_t2(:, :, mid_slice);
    slice_gt = double(volume_gt(:, :, mid_slice));

    slice_norm = (slice_t2 - min(slice_t2(:))) / (max(slice_t2(:)) - min(slice_t2(:)));
    binary_gt = slice_gt > 0;

    filtered_slice = medfilt2(slice_norm, [3 3]);
    level = graythresh(filtered_slice);
    initial_mask = imbinarize(filtered_slice, level);
    clean_mask = bwareafilt(imfill(initial_mask, 'holes'), 1);
    enhanced_slice = apply_bcet(filtered_slice, clean_mask);

    num_clusters = 4;
    [segmented_slice, candidate_tumor_mask] = apply_fcm_clustering(enhanced_slice, clean_mask, num_clusters);
    final_tumor_mask = isolate_tumor_mass(candidate_tumor_mask);
    [tumor_edges, clinical_overlay] = extract_tumor_edges(final_tumor_mask, slice_norm);

    display_preliminary_visualization(slice_norm, slice_gt, filtered_slice, ...
        clean_mask, enhanced_slice, num_clusters, segmented_slice, ...
        candidate_tumor_mask, final_tumor_mask, tumor_edges, clinical_overlay, binary_gt, plots_dir);
end

function process_optimized_patient(data_dir, patient_id, plots_dir)
% Mirrors opt_src/opt_main.m's pipeline for a single patient, saving figures only.
    image_path = fullfile(data_dir, 'imagesTr', patient_id);
    label_path = fullfile(data_dir, 'labelsTr', patient_id);

    volume4D = niftiread(niftiinfo(image_path));
    volume_gt = double(niftiread(niftiinfo(label_path)));
    volume_flair = double(volume4D(:, :, :, 1));

    brain_voxels = volume_flair(volume_flair > 0);
    hyper_threshold = quantile(brain_voxels, 0.95);
    hyperintense_area_per_slice = squeeze(sum(sum(volume_flair > hyper_threshold, 1), 2));
    [~, best_slice_idx] = max(hyperintense_area_per_slice);

    slice_flair = volume_flair(:, :, best_slice_idx);
    slice_gt = volume_gt(:, :, best_slice_idx);
    slice_norm = (slice_flair - min(slice_flair(:))) / (max(slice_flair(:)) - min(slice_flair(:)));
    binary_gt = slice_gt > 0;

    filtered_slice = medfilt2(slice_norm, [3 3]);
    level = graythresh(filtered_slice);
    initial_mask = imbinarize(filtered_slice, level);
    clean_mask = bwareafilt(imfill(initial_mask, 'holes'), 1);
    enhanced_slice = apply_bcet(filtered_slice, clean_mask);

    num_clusters = 4;
    [segmented_slice, candidate_tumor_mask] = apply_fcm_clustering(enhanced_slice, clean_mask, num_clusters);
    final_tumor_mask = isolate_tumor_mass(candidate_tumor_mask);
    [tumor_edges, clinical_overlay] = extract_tumor_edges(final_tumor_mask, slice_norm);

    display_preliminary_visualization(slice_norm, slice_gt, filtered_slice, ...
        clean_mask, enhanced_slice, num_clusters, segmented_slice, ...
        candidate_tumor_mask, final_tumor_mask, tumor_edges, clinical_overlay, binary_gt, plots_dir);
end
