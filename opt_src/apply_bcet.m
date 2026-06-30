function output_image = apply_bcet(input_image, mask)
    valid_pixels = input_image(mask);
    l = min(valid_pixels(:));           
    h = max(valid_pixels(:));           
    e = mean(valid_pixels(:));          
    s = mean(valid_pixels(:).^2);       
    L = 0.0; H = 1.0; E = 0.5;  
    denom = 2 * (h * (E - L) - e * (H - L) + l * (H - E));
    if denom == 0
        output_image = input_image; return;
    end
    num = (h^2) * (E - L) - s * (H - L) + (l^2) * (H - E);
    b = num / denom;
    a = (H - L) / ((h - b)^2 - (l - b)^2);
    c = L - a * (l - b)^2;
    output_image = zeros(size(input_image));
    output_image(mask) = a * (input_image(mask) - b).^2 + c;
    output_image(output_image < 0) = 0;
    output_image(output_image > 1) = 1;
end