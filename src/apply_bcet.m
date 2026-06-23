function output_image = apply_bcet(input_image, mask)
% APPLY_BCET Applies the Balance Contrast Enhancement Technique.
%   output_image = apply_bcet(input_image, mask) computes the parabolic
%   transformation to stretch the contrast of the input image, calculating
%   the statistical parameters only on the region defined by 'mask'.

% Extract valid pixels to compute statistics safely
valid_pixels = input_image(mask);

% Input statistics (from valid region only)
l = min(valid_pixels(:));           % Minimum input value
h = max(valid_pixels(:));           % Maximum input value
e = mean(valid_pixels(:));          % Mean input value
s = mean(valid_pixels(:).^2);       % Mean square input value

% Desired output statistics (normalized space [0, 1])
L = 0.0;  % Desired minimum
H = 1.0;  % Desired maximum
E = 0.5;  % Desired mean (balanced mid-gray)

% Compute the parabolic coefficients: y = a(x - b)^2 + c
% Denominator for parameter 'b'
denom = 2 * (h * (E - L) - e * (H - L) + l * (H - E));

% Prevent division by zero in homogeneous images
if denom == 0
    output_image = input_image;
    warning('BCET denominator is zero. Returning original image.');
    return;
end

% Numerator for parameter 'b'
num = (h^2) * (E - L) - s * (H - L) + (l^2) * (H - E);

% Calculate coefficients
b = num / denom;
a = (H - L) / ((h - b)^2 - (l - b)^2);
c = L - a * (l - b)^2;

% Initialize output image
output_image = zeros(size(input_image));

% Apply transformation only to the masked pixels
% Formula: y = a * (x - b)^2 + c
output_image(mask) = a * (input_image(mask) - b).^2 + c;

% Clamp values to ensure they stay strictly within [0, 1] range
output_image(output_image < 0) = 0;
output_image(output_image > 1) = 1;

end