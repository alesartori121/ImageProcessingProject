function fig = render_brain_tumor_3d(brain_mask_3D, tumor_mask_3D, voxel_spacing, patient_name, visible_state)
    if nargin < 5 || isempty(visible_state)
        visible_state = 'on';
    end

    smoothed_brain = smooth3(double(brain_mask_3D), 'box', 5);
    fv_brain = isosurface(smoothed_brain, 0.5);
    fv_brain.vertices = fv_brain.vertices(:, [2 1 3]) .* voxel_spacing;

    fig = figure('Name', ['3-D Tumor Reconstruction: ', patient_name], ...
                 'Position', [100, 100, 900, 750], 'Visible', visible_state);
    hold on;
    p_brain = patch(fv_brain, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.12);

    legend_handles = p_brain;
    legend_labels = {'Brain surface'};
    if nnz(tumor_mask_3D) > 0
        smoothed_tumor = smooth3(double(tumor_mask_3D), 'gaussian', 11, 2.5);
        fv_tumor = isosurface(smoothed_tumor, 0.5);
        fv_tumor.vertices = fv_tumor.vertices(:, [2 1 3]) .* voxel_spacing;
        p_tumor = patch(fv_tumor, 'FaceColor', [0.85 0.1 0.1], 'EdgeColor', 'none', 'FaceAlpha', 1.0);
        legend_handles(end+1) = p_tumor;
        legend_labels{end+1} = 'Reconstructed tumor mass';
    end

    daspect([1 1 1]);
    view(3);
    axis tight; grid on; box on;
    xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
    title(['3-D Reconstruction: Brain Shell & Tumor Mass (', patient_name, ')']);
    camlight('headlight'); camlight('left'); lighting gouraud;
    legend(legend_handles, legend_labels, 'Location', 'northeastoutside');
    hold off;
end
