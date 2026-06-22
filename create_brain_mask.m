function brain_mask = create_brain_mask(img_double)
% CREATE_BRAIN_MASK - Generates a binary mask of the brain
%   brain_mask = create_brain_mask(img_double)
%
%   Uses Otsu's method and morphological operations to isolate the brain
%   tissue from the background, improving subsequent analysis.

% Calculate a low threshold to separate the head from the black background
threshold = graythresh(img_double) * 0.2; 
mask = imbinarize(img_double, threshold);

% Morphological cleaning of the mask
mask = bwareafilt(mask, 1); % Keep ONLY the largest object (the brain)
brain_mask = imfill(mask, 'holes'); % Fill any internal gaps
end