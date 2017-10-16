function [xCenter, yCenter, mask] = calculateCenter(image, maskFactor, maskMinimum, fraction)

if maskMinimum == 0
    minIntensity = min(image(:));
else
    minIntensity = prctile(image(:), maskMinimum);
end;
meanIntensity = mean(image(image>minIntensity));

level = minIntensity + (meanIntensity - minIntensity) * maskFactor;

mask = image > level;
mask = bwareaopen(mask, round(numel(image) * fraction));

[xCoordinates, yCoordinates] = ind2sub(size(image), find(mask));
xCenter = mean(xCoordinates);
yCenter = mean(yCoordinates);

end