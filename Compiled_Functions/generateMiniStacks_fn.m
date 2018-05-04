function generateMiniStacks_fn(filename)
outputTypes = {'klb', 'jp2', 'tif'};

input_parameters = loadjson(fileread(filename));
input_parameters = convert_limits_to_values(input_parameters, {'timepoints'});
if ~isfield(input_parameters, 'verbose'), input_parameters.verbose = true; end
if input_parameters.verbose
    disp(['Using input file: ' filename]);
    disp(input_parameters)
end

splitFilename = strsplit(filename,{'_','.'});
stepName=splitFilename{end-1};

switch stepName
    case "clusterPT"
        if input_parameters.inputType == 4
            if ~isempty(input_parameters.outputLabel)
                input_parameters.outputFolder = [input_parameters.inputFolder '.corrected.' input_parameters.outputLabel];
            else
                input_parameters.outputFolder = [input_parameters.inputFolder '.corrected'];
            end
        else
            if ~isempty(input_parameters.outputLabel)
                input_parameters.outputFolder = [input_parameters.inputFolder '.corrected.' input_parameters.outputLabel filesep 'SPM' num2str(input_parameters.specimen, '%.2d')];
            else
                input_parameters.outputFolder = [input_parameters.inputFolder '.corrected' filesep 'SPM' num2str(input_parameters.specimen, '%.2d')];
            end
        end
        genericOutputLocation = [input_parameters.outputFolder filesep 'TM??????'];
    case "clusterMF"
        genericOutputLocation = [input_parameters.outputString '.TM??????_multiFused' input_parameters.outputID];
    case "clusterTF"
        genericOutputLocation = [input_parameters.outputString '.TM??????_timeFused' input_parameters.outputID];
    case "clusterCS"
        genericOutputLocation = [input_parameters.outputRoot filesep '' input_parameters.headerPattern];
    case "clusterFR"
        genericOutputLocation = [input_parameters.outputDir input_parameters.header '.TM??????' input_parameters.footer];
end

genericOutputLocation=strrep(genericOutputLocation, 'TM??????', sprintf('TM%06d', input_parameters.timepoints(1)));
if genericOutputLocation(end)==filesep, genericOutputLocation(end)=[]; end
relevantFileNames = dir([genericOutputLocation filesep '*.' outputTypes{input_parameters.outputType+1}]);
numTimepoints = numel(input_parameters.timepoints);
initialTimepoint = input_parameters.timepoints(1);
outputDirectory = [fileparts(genericOutputLocation) '/miniStacks/'];
mkdir(outputDirectory);
for fileIndex=1:numel(relevantFileNames)
    maxSlices=0; maxResizedSlices=0;
    currentRelevantFile = relevantFileNames(fileIndex);
    templateFileName = [currentRelevantFile.folder filesep  currentRelevantFile.name];
    miniStackOutputFileName = strrep(currentRelevantFile.name, sprintf('_TM%06d_', input_parameters.timepoints(1)),'_');
    miniStackOutput = [outputDirectory '/miniStack_' miniStackOutputFileName(1:end-4)...
        '.' sprintf('TM%06d', initialTimepoint) '-' sprintf('TM%06d', input_parameters.timepoints(end)) '.tif'];
    allImages=cell(numTimepoints,1);
    for index=1:numTimepoints
        timepoint=input_parameters.timepoints(index);
        im = readImage(strrep(templateFileName, sprintf('TM%06d',initialTimepoint), sprintf('TM%06d',timepoint)));
        imageSize=size(im);
        minimumDimension=min(imageSize(1:2));
        if isfield(input_parameters,'rescaleFactor')
            rescaleFactor = input_parameters.rescaleFactor;
        else
            rescaleFactor = min(200/minimumDimension, 1);
        end
        im=imresize3(im,rescaleFactor);
        resizedImageSize=size(im);
        if ismatrix(im)
            currentIm.slices = 1;
            currentIm.resizedSlices = 1;
        else
            currentIm.slices = imageSize(3);
            currentIm.resizedSlices = resizedImageSize(3);
        end

        currentIm.resizedIm = im;
        allImages{index} = currentIm;
        maxSlices = max(maxSlices, currentIm.slices);
        maxResizedSlices = max(maxResizedSlices, currentIm.resizedSlices);
    end
    
    delete(miniStackOutput);
    for count=1:numTimepoints
        currentIm = allImages{count};
        if maxResizedSlices==1
            imwrite(uint16(currentIm.resizedIm),miniStackOutput,'writemode','append');
        else
            sliceTranslation = ((maxSlices-currentIm.slices)/2)*(maxResizedSlices/maxSlices);
            imageToTranslate = zeros(size(im,1), size(im,2), maxResizedSlices);
            imageToTranslate(:,:, 1:currentIm.resizedSlices) = currentIm.resizedIm;
            translatedImage = imtranslate(imageToTranslate, [0,0,sliceTranslation]);
            for slice=1:maxResizedSlices
                imwrite(uint16(translatedImage(:,:,slice)),miniStackOutput,'writemode','append');
            end
        end
    end
    if maxResizedSlices==1
        imageDescription = sprintf('ImageJ=1.43d\nimages=%d\nframes=%d',numTimepoints, numTimepoints);
    else
        imageDescription = sprintf('ImageJ=1.43d\nimages=%d\nslices=%d\nframes=%d\nhyperstack=true\nmode=composite',maxResizedSlices*numTimepoints, maxResizedSlices,numTimepoints);
    end
    t = Tiff(miniStackOutput,'r+');
    for count=1:numTimepoints*maxResizedSlices
        setDirectory(t,count)
        setTag(t,Tiff.TagID.ImageDescription, imageDescription);
        rewriteDirectory(t);
    end
    close(t)
end
end