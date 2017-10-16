function [xOffset, yOffset] = correlatePhases(image1, image2)

if size(image1) ~= size(image2)
    error('Input images must have the same size.');
end;

filter1 = hamming(size(image1, 1)); % set up window in x-direction
filter2 = hamming(size(image1, 2)); % set up window in y-direction
imageFilter = (filter1 * filter2'); % generate 2D window function matched to image size

F1 = fft2(double(image1) .* imageFilter);
F2 = fft2(double(image2) .* imageFilter);

pdm = exp(1i * (angle(F1) - angle(F2)));
pcf = real(ifft2(pdm));
pcf = medfilt2(pcf);

[xOffset, yOffset] = find(pcf == max(pcf(:)), 1);
xOffset = xOffset - 1;
yOffset = yOffset - 1;

if xOffset > (size(image1, 1) / 2)
    xOffset = xOffset - size(image1, 1);
end;

if yOffset > (size(image1, 2) / 2)
    yOffset = yOffset - size(image1, 2);
end;

end