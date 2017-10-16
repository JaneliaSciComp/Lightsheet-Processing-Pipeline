timePoints    = [1 125 250];
referenceType = 'fusedStack.b123';
referenceTau  = 12;
dataTypes     = {...
    'CBFLarge.b300';...
    'CBFSmall.b300';...
    'MVDLarge.b300';...
    'MVDSmall.b300'};
taus          = {...
    [50 100 150 200 300 500];...
    [50 100 150 200 300 500];...
    [50 100 150 200 300 500];...
    [50 100 150 200 300 500]};
rois          = {...
    [149, 149+1208, 517, 517+648, 389, 1101];...
    [76, 76+1494, 108, 108+1096, 278, 1191];...
    [7, 7+1784, 26, 26+1315, 69, 1418]};
threshold     = 250;
background    = 100;

for t = 1:numel(timePoints)
    imageStack = readImage(['TM' num2str(timePoints(t), '%.6d') '.' referenceType(1:(end - 4)) 'klb']);
    imageStack = imageStack(rois{t}(1):rois{t}(2), rois{t}(3):rois{t}(4), rois{t}(5):rois{t}(6));
    if background > 0
        imageStack = imageStack - background;
    end;
    imageStack = uint8(imageStack .* (255/threshold));
    
    referenceStack = readImage(['TM' num2str(timePoints(t), '%.6d') '.' referenceType '.bin_tau' num2str(referenceTau) '.klb']);
    referenceStack = referenceStack(rois{t}(1):rois{t}(2), rois{t}(3):rois{t}(4), rois{t}(5):rois{t}(6));
    referenceEdgeStack = false(size(referenceStack));
    for z = 1:size(referenceStack, 3)
        referenceEdgeStack(:, :, z) = imdilate(referenceStack(:, :, z), ones(3, 3)) > imerode(referenceStack(:, :, z), ones(3, 3));
        % referenceEdgeStack(:, :, z) = edge(referenceStack(:, :, z), 'log', 0);
    end;
    referenceEdgeStack = uint8(referenceEdgeStack) .* 255;
    
    for i = 1:numel(dataTypes)
        for j = 1:numel(taus{i})
            segmentationStack = readImage(['TM' num2str(timePoints(t), '%.6d') '.' dataTypes{i} '.bin_tau' num2str(taus{i}(j)) '.klb']);
            segmentationStack = segmentationStack(rois{t}(1):rois{t}(2), rois{t}(3):rois{t}(4), rois{t}(5):rois{t}(6));
            edgeStack = false(size(segmentationStack));
            for z = 1:size(segmentationStack, 3)
                edgeStack(:, :, z) = imdilate(segmentationStack(:, :, z), ones(3, 3)) > imerode(segmentationStack(:, :, z), ones(3, 3));
                % edgeStack(:, :, z) = edge(segmentationStack(:, :, z), 'log', 0);
            end;
            edgeStack = uint8(edgeStack) .* 255;
            
            folderName = ['Comparison.TM' num2str(timePoints(t), '%.6d') '.' dataTypes{i} '.bin_tau' num2str(taus{i}(j)) '.overlay'];
            mkdir(folderName);
            for z = 1:size(imageStack, 3)
                imwrite(cat(3, edgeStack(:, :, z), referenceEdgeStack(:, :, z), imageStack(:, :, z)), [folderName '\' num2str(z, '%.4d') '.tif'], 'Compression', 'none');
            end;
        end;
    end;
end;