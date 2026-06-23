function [edge_map, overlay_img] = extract_tumor_edges(tumor_mask, original_slice)
% EXTRACT_TUMOR_EDGES Applies Canny detector and generates clinical overlay.
    edge_map = edge(tumor_mask, 'canny');
    
    % Create RGB overlay
    overlay_img = repmat(original_slice, [1, 1, 3]);
    r = overlay_img(:, :, 1); g = overlay_img(:, :, 2); b = overlay_img(:, :, 3);
    
    r(edge_map) = 1.0; g(edge_map) = 0.0; b(edge_map) = 0.0; % Red edges
    overlay_img(:, :, 1) = r; overlay_img(:, :, 2) = g; overlay_img(:, :, 3) = b;
end