function filterResults(timepoints, inputDir, outputDir, header, footer, stackLabel, fusionFlag, specimen, cameras, channels,subOffset, ...
                        removeDirt, filterMode, rangeArray, splitting, scaling, preMedian, postMedian, inputType, outputType, ...
                        subProject, saveRawMax, saveStacks, jobMemory, t, memoryEstimate)

% -----------------------------------------------------------------------------------------------
% | Adaptive image background correction                                                        |
% |                                                                                             |
% | Code by Philipp J. Keller, HHMI/Janelia Research Campus, 2011-2017                          |
% | Email: kellerp@janelia.hhmi.org                                                             |
% |                                                                                             |
% | Utilizes optimization modules and functions by Fernando Amat, HHMI/Janelia Research Campus: |
% | imgaussianAnisotropy.m (included below)                                                     |
% | readKLBstack.mexw64                                                                         |
% | writeKLBstack.mexw64                                                                        |
% -----------------------------------------------------------------------------------------------

timepoint = timepoints(t);

version = 1.10;

configuration = cell(25, 1);

configuration{1}  = version;      configuration{2}  = timepoint;   configuration{3}  = inputDir;    configuration{4}  = outputDir;
configuration{5}  = header;       configuration{6}  = footer;      configuration{7}  = stackLabel;
configuration{8}  = fusionFlag;
configuration{9}  = specimen;     configuration{10} = cameras;     configuration{11} = channels;
configuration{12} = subOffset;
configuration{13} = removeDirt;
configuration{14} = filterMode;   configuration{15} = rangeArray;  configuration{16} = splitting;   configuration{17} = scaling;
configuration{18} = preMedian;    configuration{19} = postMedian;  configuration{20} = inputType;   configuration{21} = outputType;
configuration{22} = subProject;   configuration{23} = saveRawMax;  configuration{24} = saveStacks;
configuration{25} = [jobMemory(1) memoryEstimate];         

if fusionFlag == 0
    configurationString = ['_CM' num2str(cameras(1), '%.2d') '_CHN' num2str(channels(1), '%.2d')]; % single channel
elseif length(cameras) == 2 && length(channels) == 2
    configurationString = ['_CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(channels(1), '%.2d') '_CHN' num2str(channels(2), '%.2d')]; % 4-view fusion
elseif length(cameras) == 1 && length(channels) == 2
    configurationString = ['_CM' num2str(cameras(1), '%.2d') '_CHN' num2str(channels(1), '%.2d') '_CHN' num2str(channels(2), '%.2d')]; % 2-view channel fusion
elseif length(cameras) == 2 && length(channels) == 1
    configurationString = ['_CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(channels(1), '%.2d')]; % 2-view camera fusion
else
    error('Error: Provide either 2 channels and 2 cameras (4-view fusion) or 2 channels and 1 camera (2-view channel fusion) or 1 channel and 2 cameras (2-view camera fusion) or unfused data (fusionFlag == 0)');
end;

switch inputType
    case 0
        inputExtension = '.klb';
    case 1
        inputExtension = '.jp2';
    case 2
        inputExtension = '.tif';
end;

switch outputType
    case 0
        outputExtension = '.klb';
    case 1
        outputExtension = '.jp2';
    case 2
        outputExtension = '.tif';
end;

disp(['processing time point ' num2str(timepoint, '%.6d')]);

fileNameHeader = ['SPM' num2str(specimen, '%.2d') '_TM' num2str(timepoint, '%.6d') configurationString];
if fusionFlag == 0
    stackName = [inputDir filesep 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timepoint, '%.6d') filesep fileNameHeader inputExtension];
else
    stackName = [inputDir filesep header '.TM' num2str(timepoint, '%.6d') footer filesep '' fileNameHeader '.fusedStack' stackLabel inputExtension];
end

stack = readImage(stackName);
xSize = size(stack, 1);
ySize = size(stack, 2);
zSize = size(stack, 3);

if subOffset ~= 0
    stack = stack - subOffset;
end;

if fusionFlag == 0
    outputPath = [outputDir 'SPM' num2str(specimen, '%.2d') filesep 'TM' num2str(timepoint, '%.6d')];
    intermediateString = '.';
    intermediateStackString = '';
else
    outputPath = [outputDir filesep header '.TM' num2str(timepoint, '%.6d') footer];
    intermediateString = intermediateString '';
    intermediateStackString = intermediateStackString '';
end
mkdir(outputPath);

save([outputPath filesep '' fileNameHeader '.configuration.mat'], 'configuration');

if subProject
    center = round(zSize / 2);
    
    firstHalfProjected  = max(stack(:, :, 1:center), [], 3);
    secondHalfProjected = max(stack(:, :, (center + 1):end), [], 3);
    
    firstHalfProjectedName  = [outputPath filesep '' fileNameHeader intermediateString 'xyFProjection' outputExtension];
    secondHalfProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'fusedStack_xyBProjection' outputExtension];
    
    writeImage(firstHalfProjected, firstHalfProjectedName);
    writeImage(secondHalfProjected, secondHalfProjectedName);
end;

if saveRawMax
    xyProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'xyProjection' outputExtension];
    xzProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'xzProjection' outputExtension];
    yzProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'yzProjection' outputExtension];
    
    writeImage(max(stack, [], 3), xyProjectedName);
    writeImage(squeeze(max(stack, [], 2)), xzProjectedName);
    writeImage(squeeze(max(stack, [], 1)), yzProjectedName);
end;

if removeDirt(1)
    stackMask = stack > removeDirt(2);
    objectStats = regionprops(stackMask, 'Area', 'PixelIdxList');
    largestBlobIndex = 1;
    for i = 1:numel(objectStats)
        if objectStats(i).Area > objectStats(largestBlobIndex).Area
            largestBlobIndex = i;
        end;
    end;
    stackMask = true(size(stackMask));
    stackMask(objectStats(largestBlobIndex).PixelIdxList) = 0;
    stack(stackMask) = 0;
    
    if saveStacks
        cleanedStackName = [outputPath filesep '' fileNameHeader intermediateStackString '.cleaned' outputExtension];
        writeImage(stack, cleanedStackName);
    end;
    
    if subProject
        firstHalfProjected  = max(stack(:, :, 1:center), [], 3);
        secondHalfProjected = max(stack(:, :, (center + 1):end), [], 3);
        
        firstHalfProjectedName  = [outputPath filesep '' fileNameHeader intermediateString 'xyFProjection.cleaned' outputExtension];
        secondHalfProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'xyBProjection.cleaned' outputExtension];
        
        writeImage(firstHalfProjected, firstHalfProjectedName);
        writeImage(secondHalfProjected, secondHalfProjectedName);
    end;
    
    xyProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'xyProjection.cleaned' outputExtension];
    xzProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'xzProjection.cleaned' outputExtension];
    yzProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'yzProjection.cleaned' outputExtension];
    
    writeImage(max(stack, [], 3), xyProjectedName);
    writeImage(squeeze(max(stack, [], 2)), xzProjectedName);
    writeImage(squeeze(max(stack, [], 1)), yzProjectedName);
end;

if preMedian(1)
    stackCopy = stack;
    for z = 1:size(stackCopy, 3)
        stackCopy(:, :, z) = medfilt2(stackCopy(:, :, z), [preMedian(2) preMedian(2)]);
    end;
    
    if saveStacks
        if removeDirt(1)
            filteredStackName = [outputPath filesep '' fileNameHeader intermediateStackString '.cleaned_median' outputExtension];
        else
            filteredStackName = [outputPath filesep '' fileNameHeader intermediateStackString '.median' outputExtension];
        end;
        writeImage(stackCopy, filteredStackName);
    end;
    
    if subProject
        firstHalfFilteredProjected  = max(stackCopy(:, :, 1:center), [], 3);
        secondHalfFilteredProjected = max(stackCopy(:, :, (center + 1):end), [], 3);
        
        if removeDirt(1)
            firstHalfFilteredProjectedName  = [outputPath filesep '' fileNameHeader intermediateString 'xyFProjection.cleaned_median' outputExtension];
            secondHalfFilteredProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'xyBProjection.cleaned_median' outputExtension];
        else
            firstHalfFilteredProjectedName  = [outputPath filesep '' fileNameHeader intermediateString 'xyFProjection.median' outputExtension];
            secondHalfFilteredProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'xyBProjection.median' outputExtension];
        end;
        
        writeImage(firstHalfFilteredProjected, firstHalfFilteredProjectedName);
        writeImage(secondHalfFilteredProjected, secondHalfFilteredProjectedName);
    end;
    
    if removeDirt(1)
        xyFilteredProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'xyProjection.cleaned_median' outputExtension];
        xzFilteredProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'xzProjection.cleaned_median' outputExtension];
        yzFilteredProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'yzProjection.cleaned_median' outputExtension];
    else
        xyFilteredProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'xyProjection.median' outputExtension];
        xzFilteredProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'xzProjection.median' outputExtension];
        yzFilteredProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'yzProjection.median' outputExtension];
    end;
    
    writeImage(max(stackCopy, [], 3), xyFilteredProjectedName);
    writeImage(squeeze(max(stackCopy, [], 2)), xzFilteredProjectedName);
    writeImage(squeeze(max(stackCopy, [], 1)), yzFilteredProjectedName);
    
    clear stackCopy;
end;

if ~isempty(filterMode)
    for r = rangeArray
        kernelSizeArray  = [r r max(1, r / scaling)];
        kernelSigmaArray = [r r max(1, r / scaling)];
        
        switch filterMode
            case 0
                filteredStack = stack - medfilt3(stack, r);
            case 1
                filteredStack = stack - uint16(smooth3(stack, 'box', r));
            case 2
                if splitting == 0 || splitting == 1
                    filteredStack = stack - uint16(imgaussianAnisotropy(double(stack), kernelSizeArray, kernelSigmaArray));
                else
                    gaussStack = zeros(xSize, ySize, zSize, 'uint16');
                    splittingMargin = r + 1;
                    
                    for i = 1:splitting
                        xSlabStart = max(1, round((i - 1) * xSize / splitting + 1 - splittingMargin));
                        xSlabStop = min(xSize, round(i * xSize / splitting + splittingMargin));
                        convolvedSlab = uint16(imgaussianAnisotropy(double(stack(xSlabStart:xSlabStop, :, :)), kernelSizeArray, kernelSigmaArray));
                        if i == 1
                            gaussStack(1:(xSlabStop - splittingMargin), :, :) = convolvedSlab(1:(end - splittingMargin), :, :);
                        elseif i == splitting
                            gaussStack((xSlabStart + splittingMargin):end, :, :) = convolvedSlab((1 + splittingMargin):end, :, :);
                        else % i > 1 && i < splitting
                            gaussStack((xSlabStart + splittingMargin):(xSlabStop - splittingMargin), :, :) = convolvedSlab((1 + splittingMargin):(end - splittingMargin), :, :);
                        end;
                        clear convolvedSlab;
                    end;
                    
                    filteredStack = stack - gaussStack;
                    clear gaussStack;
                end;
        end;
        
        if saveStacks
            filteredStackName = [outputPath filesep '' fileNameHeader intermediateStackString '.filtered_' num2str(r) outputExtension];
            writeImage(filteredStack, filteredStackName);
        end;
        
        if subProject
            firstHalfFilteredProjected  = max(filteredStack(:, :, 1:center), [], 3);
            secondHalfFilteredProjected = max(filteredStack(:, :, (center + 1):end), [], 3);
            
            firstHalfFilteredProjectedName  = [outputPath filesep '' fileNameHeader intermediateString 'xyFProjection.filtered_' num2str(r) outputExtension];
            secondHalfFilteredProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'xyBProjection.filtered_' num2str(r) outputExtension];
            
            writeImage(firstHalfFilteredProjected, firstHalfFilteredProjectedName);
            writeImage(secondHalfFilteredProjected, secondHalfFilteredProjectedName);
        end;
        
        xyFilteredProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'xyProjection.filtered_' num2str(r) outputExtension];
        xzFilteredProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'xzProjection.filtered_' num2str(r) outputExtension];
        yzFilteredProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'yzProjection.filtered_' num2str(r) outputExtension];
        
        writeImage(max(filteredStack, [], 3), xyFilteredProjectedName);
        writeImage(squeeze(max(filteredStack, [], 2)), xzFilteredProjectedName);
        writeImage(squeeze(max(filteredStack, [], 1)), yzFilteredProjectedName);
        
        if postMedian(1)
            for z = 1:size(filteredStack, 3)
                filteredStack(:, :, z) = medfilt2(filteredStack(:, :, z), [postMedian(2) postMedian(2)]);
            end;
            
            if saveStacks
                filteredStackName = [outputPath filesep '' fileNameHeader intermediateStackString '.filtered_' num2str(r) '_median' outputExtension];
                writeImage(filteredStack, filteredStackName);
            end;
            
            if subProject
                firstHalfFilteredProjected  = max(filteredStack(:, :, 1:center), [], 3);
                secondHalfFilteredProjected = max(filteredStack(:, :, (center + 1):end), [], 3);
                
                firstHalfFilteredProjectedName  = [outputPath filesep '' fileNameHeader intermediateString 'xyFProjection.filtered_' num2str(r) '_median' outputExtension];
                secondHalfFilteredProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'xyBProjection.filtered_' num2str(r) '_median' outputExtension];
                
                writeImage(firstHalfFilteredProjected, firstHalfFilteredProjectedName);
                writeImage(secondHalfFilteredProjected, secondHalfFilteredProjectedName);
            end;
            
            xyFilteredProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'xyProjection.filtered_' num2str(r) '_median' outputExtension];
            xzFilteredProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'xzProjection.filtered_' num2str(r) '_median' outputExtension];
            yzFilteredProjectedName = [outputPath filesep '' fileNameHeader intermediateString 'yzProjection.filtered_' num2str(r) '_median' outputExtension];
            
            writeImage(max(filteredStack, [], 3), xyFilteredProjectedName);
            writeImage(squeeze(max(filteredStack, [], 2)), xzFilteredProjectedName);
            writeImage(squeeze(max(filteredStack, [], 1)), yzFilteredProjectedName);
        end;
        
        clear filteredStack;
    end;
end;

clear stack;

end

function I = imgaussianAnisotropy(I, sigma, siz)

% --------------------------------------------------------------------------------------------------------------
% | Anisotropic Gaussian filtering                                                                             |
% | Original function (imgaussian) written by D. Kroon, University of Twente, September 2009                   |
% | Code modification to allow different sigmas by Fernando Amat, HHMI/Janelia Research Campus, September 2010 |
% --------------------------------------------------------------------------------------------------------------

% Note: If X is of type double, the MEX file is faster. For images of type single or int, this m-file is faster.

if(~exist('siz', 'var'))
    siz = sigma * 6;
end;

ndimsI = sum(size(I) > 1);
if(length(sigma) ~= ndimsI)
    error 'You must specify one sigma for each dimension of the image'
end;

% Filter each dimension with the 1D Gaussian kernels
if(ndimsI == 1)
    % Make 1D Gaussian kernel
    kk = 1;
    x = -ceil(siz(kk) / 2):ceil(siz(kk) / 2);
    H = exp(-(x .^ 2 / (2 * sigma(kk) ^ 2)));
    H = H / sum(H(:));

    I = imfilter(I, H, 'same', 'replicate');
elseif(ndimsI == 2)
    % Make 1D Gaussian kernel
    kk = 1;
    x = -ceil(siz(kk) / 2):ceil(siz(kk) / 2);
    H = exp(-(x .^ 2 / (2 * sigma(kk) ^ 2)));
    H = H / sum(H(:));
    Hx = reshape(H, [length(H) 1]);
    
    kk = 2;
    x = -ceil(siz(kk) / 2):ceil(siz(kk) / 2);
    H = exp(-(x .^ 2 / (2 * sigma(kk) ^ 2)));
    H = H / sum(H(:));
    Hy = reshape(H, [1 length(H)]);
    
    I = imfilter(imfilter(I, Hx, 'same', 'replicate'), Hy, 'same', 'replicate');
elseif(ndimsI == 3)
    if(size(I, 3) < 4) % Detect if 3D or color image
        % Make 1D Gaussian kernel
        kk = 1;
        x = -ceil(siz(kk) / 2):ceil(siz(kk) / 2);
        H = exp(-(x .^ 2 / (2 * sigma(kk) ^ 2)));
        H = H / sum(H(:));
        Hx = reshape(H, [length(H) 1]);
        
        kk = 2;
        x = -ceil(siz(kk) / 2):ceil(siz(kk) / 2);
        H = exp(-(x .^ 2 / (2 * sigma(kk) ^ 2)));
        H = H / sum(H(:));
        Hy = reshape(H, [1 length(H)]);
        for k = 1:size(I, 3)
            I(:, :, k) = imfilter(imfilter(I(:, :, k), Hx, 'same', 'replicate'), Hy, 'same', 'replicate');
        end;
    else
        % Make 1D Gaussian kernel
        kk = 1;
        x = -ceil(siz(kk) / 2):ceil(siz(kk) / 2);
        H = exp(-(x .^ 2 / (2 * sigma(kk) ^ 2)));
        H = H / sum(H(:));
        Hx = reshape(H, [length(H) 1 1]);
        
        kk = 2;
        x = -ceil(siz(kk) / 2):ceil(siz(kk) / 2);
        H = exp(-(x .^ 2 / (2 * sigma(kk) ^ 2)));
        H = H / sum(H(:));
        Hy = reshape(H, [1 length(H) 1]);
        
        kk = 3;
        x = -ceil(siz(kk) / 2):ceil(siz(kk) / 2);
        H = exp(-(x .^ 2 / (2 * sigma(kk) ^ 2)));
        H = H / sum(H(:));
        Hz = reshape(H, [1 1 length(H)]);
        
        I = imfilter(imfilter(imfilter(I, Hx, 'same', 'replicate'), Hy, 'same', 'replicate'), Hz, 'same', 'replicate');
    end;
else
    error('imgaussian:input', 'unsupported input dimension');
end;

end