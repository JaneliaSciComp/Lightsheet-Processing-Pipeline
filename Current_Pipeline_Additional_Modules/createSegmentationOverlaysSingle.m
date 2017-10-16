imageStackName        = 'TM000250.CBFLarge.klb';
segmentationStackName = 'TM000250.CBFLarge.b300_median2.bin_tau100.klb';
roi                   = [7, 7+1784, 26, 26+1315, 69, 1418];
background            = 0;
threshold             = 3000;

imageStack = readImage(imageStackName);
imageStack = imageStack(roi(1):roi(2), roi(3):roi(4), roi(5):roi(6));
if background > 0
    imageStack = imageStack - background;
end;
imageStack = uint8(imageStack .* (255/threshold));

segmentationStack = readImage(segmentationStackName);
segmentationStack = segmentationStack(roi(1):roi(2), roi(3):roi(4), roi(5):roi(6));
edgeStack = false(size(segmentationStack));
for z = 1:size(segmentationStack, 3)
    edgeStack(:, :, z) = imdilate(segmentationStack(:, :, z), ones(3, 3)) > imerode(segmentationStack(:, :, z), ones(3, 3));
end;
edgeStack = uint8(edgeStack) .* 255;

emptyFrame = zeros(size(edgeStack, 1), size(edgeStack, 2), 'uint8');
folderName = [segmentationStackName(1:(end - 3)) 'overlay'];
mkdir(folderName);
for z = 1:size(imageStack, 3)
    imwrite(cat(3, edgeStack(:, :, z), imageStack(:, :, z), emptyFrame), [folderName '\' num2str(z, '%.4d') '.tif'], 'Compression', 'none');
end;