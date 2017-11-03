timePoints = [1 125 250];
dataTypes  = {...
    'CBFLarge.b300';...
    'CBFSmall.b300';...
    'fusedStack.b123';...
    'MVDLarge.b300';...
    'MVDSmall.b300'};
taus       = {...
    [50 100 150 200 300 500];...
    [50 100 150 200 300 500];...
    [10 12 15];...
    [50 100 150 200 300 500];...
    [50 100 150 200 300 500]};
thresholds = [3000, 3000, 250, 3000, 3000];
rois       = {...
    [149, 149+1208, 517, 517+648, 389, 1101];...
    [76, 76+1494, 108, 108+1096, 278, 1191];...
    [7, 7+1784, 26, 26+1315, 69, 1418]};
background = [0, 0, 100, 0, 0];

for t = 1:numel(timePoints)
    for i = 1:numel(dataTypes)
        imageStack = readImage(['TM' num2str(timePoints(t), '%.6d') '.' dataTypes{i}(1:(end - 4)) 'klb']);
        imageStack = imageStack(rois{t}(1):rois{t}(2), rois{t}(3):rois{t}(4), rois{t}(5):rois{t}(6));
        if background(i) > 0
            imageStack = imageStack - background(i);
        end;
        imageStack = uint8(imageStack .* (255/thresholds(i)));
        for j = 1:numel(taus{i})
            segmentationStack = readImage(['TM' num2str(timePoints(t), '%.6d') '.' dataTypes{i} '.bin_tau' num2str(taus{i}(j)) '.klb']);
            segmentationStack = segmentationStack(rois{t}(1):rois{t}(2), rois{t}(3):rois{t}(4), rois{t}(5):rois{t}(6));
            edgeStack = false(size(segmentationStack));
            for z = 1:size(segmentationStack, 3)
                edgeStack(:, :, z) = imdilate(segmentationStack(:, :, z), ones(3, 3)) > imerode(segmentationStack(:, :, z), ones(3, 3));
                % edgeStack(:, :, z) = edge(segmentationStack(:, :, z), 'log', 0);
            end;
            edgeStack = uint8(edgeStack) .* 255;
            
            emptyFrame = zeros(size(edgeStack, 1), size(edgeStack, 2), 'uint8');
            folderName = ['TM' num2str(timePoints(t), '%.6d') '.' dataTypes{i} '.bin_tau' num2str(taus{i}(j)) '.overlay'];
            mkdir(folderName);
            for z = 1:size(imageStack, 3)
                imwrite(cat(3, edgeStack(:, :, z), imageStack(:, :, z), emptyFrame), [folderName filesep '' num2str(z, '%.4d') '.tif'], 'Compression', 'none');
            end;
        end;
    end;
end;