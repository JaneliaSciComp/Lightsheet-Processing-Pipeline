files = dir('*.klb');
currentHeader = readKLBheader(...
    [files(end).name]);
aSize = currentHeader.xyzct(1);
bSize = currentHeader.xyzct(2);
projections = zeros(aSize, bSize, numel(files), 'uint16');

for t = 1:numel(files)
    currentProjection = readImage(...
        [files(t).name]);
    startIndex = floor((bSize - size(currentProjection, 2)) / 2) + 1;
    stopIndex = floor((bSize - size(currentProjection, 2)) / 2) + size(currentProjection, 2);
    projections(:, startIndex:stopIndex, t) = currentProjection;
end;
outputFilename = 'projections.klb';
writeImage(projections, outputFilename);