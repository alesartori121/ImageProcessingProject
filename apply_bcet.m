function img_bcet = apply_bcet(img_input)
% APPLY_BCET - Balance Contrast Enhancement Technique
%   img_bcet = apply_bcet(img_input)
%
%   Applies the Balance Contrast Enhancement Technique (BCET) to an image.
%   This algorithm improves the image contrast by using a parabolic 
%   transformation while maintaining the fundamental shape of the histogram.
%
%   Reference: Zotin et al., "Edge detection in MRI brain tumor images..."

    % 1. Extract input image statistics (L, H, E)
    % We use (:) to consider the 2D matrix as a single 1D vector
    L = min(img_input(:)); % Minimum intensity
    H = max(img_input(:)); % Maximum intensity
    E = mean(img_input(:)); % Mean intensity

    % 2. Define adaptive target statistics
    % Instead of forcing absolute values that cause saturation,
    % we adapt the targets to the original dynamics of the image.
    l = 0.0; % Target minimum (dark background)
    h = min(1.0, H * 1.2); % Target maximum: slightly expand, but max 1.0
    
    % The target mean (e) is delicate. Forcing it to 0.5 on a very dark
    % image causes white saturation. We slightly increase the original mean.
    e = E * 1.2; % Increase mean by 20%
    
    % Safety check: target mean must be strictly between l and h
    if e >= h
        e = h - 0.1;
    elseif e <= l
        e = l + 0.1;
    end

    % 3. Calculate parabolic coefficients (a, b, c) based on BCET math
    numerator_b   = H^2 * (e - l) - E^2 * (h - l) + L^2 * (h - e);
    denominator_b = 2 * (H * (e - l) - E * (h - l) + L * (h - e));
    
    % Check to avoid division by zero
    if denominator_b == 0
        b = 0; 
    else
        b = numerator_b / denominator_b;
    end
    
    % Equations for 'a' and 'c'
    denom_a = (H - b)^2 - (L - b)^2;
    if denom_a == 0
        a = 1; % If image is uniform, no transformation
    else
        a = (h - l) / denom_a;
    end
    
    c = l - a * (L - b)^2;

    % 4. Apply the parabolic transformation to the entire image
    img_bcet = a .* (img_input - b).^2 + c;

    % 5. Ensure limits are respected (Clipping to avoid out-of-bounds)
    img_bcet(img_bcet < 0) = 0;
    img_bcet(img_bcet > 1) = 1;

end