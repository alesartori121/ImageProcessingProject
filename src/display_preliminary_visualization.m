function display_preliminary_visualization(slice_t2, slice_gt)
% DISPLAY_PRELIMINARY_VISUALIZATION Shows the T2 slice and Ground Truth.
%   display_preliminary_visualization(slice_t2, slice_gt) creates a
%   figure to visually inspect the normalized T2 MRI slice and the
%   corresponding Ground Truth labels before further processing.

% Create a figure with a specific size for better visibility
figure('Name', 'Data Understanding: T2 and Ground Truth', 'Position', [100, 100, 1000, 400]);

% Display the T2 slice (using [] to scale display based on min/max values)
subplot(1, 2, 1);
imshow(slice_t2, []);
title('Normalized T2 Slice');

% Display the Ground Truth labels
subplot(1, 2, 2);
imshow(slice_gt, []);
title('Original Ground Truth Labels');

% Force MATLAB to draw the figure immediately
drawnow;
end