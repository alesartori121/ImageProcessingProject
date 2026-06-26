function [edge_map, overlay_img] = extract_tumor_edges(tumor_mask, original_slice)
% EXTRACT_TUMOR_EDGES Traces the boundary of the binary tumor mask using a
% morphological edge detector, and generates a clinical overlay.
%   NOTE ON FIDELITY: tumor_mask is already binary, so a real Canny detector
%   has nothing to act on (Canny needs a grayscale gradient -- in this
%   course's slides it is only ever run on intensity images, e.g.
%   circuit.tif, headCT.tif, never on a pre-segmented binary mask). The
%   course's own worked exercise for tracing the boundary of an
%   already-binarized region (leaf2.m: "Extract a binary edge map using a
%   morphological edge detector", s=strel('square',3); xd=imdilate(x,s);
%   xc=xd-x;) uses dilation minus the original mask instead, so that is
%   used here.
    se = strel('square', 3);
    edge_map = imdilate(tumor_mask, se) & ~tumor_mask;
    
    % Create RGB overlay
    overlay_img = repmat(original_slice, [1, 1, 3]);
    r = overlay_img(:, :, 1); g = overlay_img(:, :, 2); b = overlay_img(:, :, 3);
    
    r(edge_map) = 1.0; g(edge_map) = 0.0; b(edge_map) = 0.0; % Red edges
    overlay_img(:, :, 1) = r; overlay_img(:, :, 2) = g; overlay_img(:, :, 3) = b;
end