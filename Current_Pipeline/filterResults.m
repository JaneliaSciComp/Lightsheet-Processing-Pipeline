function filterResults(parameterDatabase, t, memoryEstimate)

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

load(parameterDatabase);
timepoint = timepoints(t);

version = 1.08;

configuration = cell(23, 1);

configuration{1}  = version;      configuration{2}  = timepoint;   configuration{3}  = inputDir;    configuration{4}  = outputDir;
configuration{5}  = header;       configuration{6}  = footer;      configuration{7}  = stackLabel;
configuration{8}  = specimen;     configuration{9}  = cameras;     configuration{10} = channels;
configuration{11} = removeDirt;
configuration{12} = filterMode;   configuration{13} = rangeArray;  configuration{14} = splitting;   configuration{15} = scaling;
configuration{16} = preMedian;    configuration{17} = postMedian;  configuration{18} = inputType;   configuration{19} = outputType;
configuration{20} = subProject;   configuration{21} = saveRawMax;  configuration{22} = saveStacks;
configuration{23} = [jobMemory(1) memoryEstimate];        

if length(cameras) == 2 && length(channels) == 2
    configurationString = ['_CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(channels(1), '%.2d') '_CHN' num2str(channels(2), '%.2d')]; % 4-view fusion
elseif length(cameras) == 1 && length(channels) == 2
    configurationString = ['_CM' num2str(cameras(1), '%.2d') '_CHN' num2str(channels(1), '%.2d') '_CHN' num2str(channels(2), '%.2d')]; % 2-view channel fusion
elseif length(cameras) == 2 && length(channels) == 1
    configurationString = ['_CM' num2str(cameras(1), '%.2d') '_CM' num2str(cameras(2), '%.2d') '_CHN' num2str(channels(1), '%.2d')]; % 2-view camera fusion
else
    error('Error: Provide either 2 channels and 2 cameras (4-view fusion) or 2 channels and 1 camera (2-view channel fusion) or 1 channel and 2 cameras (2-view camera fusion)');
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
stackName = [inputDir filesep header '.TM' num2str(timepoint, '%.6d') footer filesep '' fileNameHeader '.fusedStack' stackLabel inputExtension];

stack = readImage(stackName);
xSize = size(stack, 1);
ySize = size(stack, 2);
zSize = size(stack, 3);

outputPath = [outputDir filesep header '.TM' num2str(timepoint, '%.6d') footer];
mkdir(outputPath);

save([outputPath filesep '' fileNameHeader '.configuration.mat'], 'configuration');

if subProject
    center = round(zSize / 2);
    
    firstHalfProjected  = max(stack(:, :, 1:center), [], 3);
    secondHalfProjected = max(stack(:, :, (center + 1):end), [], 3);
    
    firstHalfProjectedName  = [outputPath filesep '' fileNameHeader '.fusedStack_xyFProjection' outputExtension];
    secondHalfProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_xyBProjection' outputExtension];
    
    writeImage(firstHalfProjected, firstHalfProjectedName);
    writeImage(secondHalfProjected, secondHalfProjectedName);
end;

if saveRawMax
    xyProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_xyProjection' outputExtension];
    xzProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_xzProjection' outputExtension];
    yzProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_yzProjection' outputExtension];
    
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
        cleanedStackName = [outputPath filesep '' fileNameHeader '.fusedStack.cleaned' outputExtension];
        writeImage(stack, cleanedStackName);
    end;
    
    if subProject
        firstHalfProjected  = max(stack(:, :, 1:center), [], 3);
        secondHalfProjected = max(stack(:, :, (center + 1):end), [], 3);
        
        firstHalfProjectedName  = [outputPath filesep '' fileNameHeader '.fusedStack_xyFProjection.cleaned' outputExtension];
        secondHalfProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_xyBProjection.cleaned' outputExtension];
        
        writeImage(firstHalfProjected, firstHalfProjectedName);
        writeImage(secondHalfProjected, secondHalfProjectedName);
    end;
    
    xyProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_xyProjection.cleaned' outputExtension];
    xzProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_xzProjection.cleaned' outputExtension];
    yzProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_yzProjection.cleaned' outputExtension];
    
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
            filteredStackName = [outputPath filesep '' fileNameHeader '.fusedStack.cleaned_median' outputExtension];
        else
            filteredStackName = [outputPath filesep '' fileNameHeader '.fusedStack.median' outputExtension];
        end;
        writeImage(stackCopy, filteredStackName);
    end;
    
    if subProject
        firstHalfFilteredProjected  = max(stackCopy(:, :, 1:center), [], 3);
        secondHalfFilteredProjected = max(stackCopy(:, :, (center + 1):end), [], 3);
        
        if removeDirt(1)
            firstHalfFilteredProjectedName  = [outputPath filesep '' fileNameHeader '.fusedStack_xyFProjection.cleaned_median' outputExtension];
            secondHalfFilteredProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_xyBProjection.cleaned_median' outputExtension];
        else
            firstHalfFilteredProjectedName  = [outputPath filesep '' fileNameHeader '.fusedStack_xyFProjection.median' outputExtension];
            secondHalfFilteredProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_xyBProjection.median' outputExtension];
        end;
        
        writeImage(firstHalfFilteredProjected, firstHalfFilteredProjectedName);
        writeImage(secondHalfFilteredProjected, secondHalfFilteredProjectedName);
    end;
    
    if removeDirt(1)
        xyFilteredProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_xyProjection.cleaned_median' outputExtension];
        xzFilteredProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_xzProjection.cleaned_median' outputExtension];
        yzFilteredProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_yzProjection.cleaned_median' outputExtension];
    else
        xyFilteredProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_xyProjection.median' outputExtension];
        xzFilteredProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_xzProjection.median' outputExtension];
        yzFilteredProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_yzProjection.median' outputExtension];
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
            filteredStackName = [outputPath filesep '' fileNameHeader '.fusedStack.filtered_' num2str(r) outputExtension];
            writeImage(filteredStack, filteredStackName);
        end;
        
        if subProject
            firstHalfFilteredProjected  = max(filteredStack(:, :, 1:center), [], 3);
            secondHalfFilteredProjected = max(filteredStack(:, :, (center + 1):end), [], 3);
            
            firstHalfFilteredProjectedName  = [outputPath filesep '' fileNameHeader '.fusedStack_xyFProjection.filtered_' num2str(r) outputExtension];
            secondHalfFilteredProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_xyBProjection.filtered_' num2str(r) outputExtension];
            
            writeImage(firstHalfFilteredProjected, firstHalfFilteredProjectedName);
            writeImage(secondHalfFilteredProjected, secondHalfFilteredProjectedName);
        end;
        
        xyFilteredProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_xyProjection.filtered_' num2str(r) outputExtension];
        xzFilteredProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_xzProjection.filtered_' num2str(r) outputExtension];
        yzFilteredProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_yzProjection.filtered_' num2str(r) outputExtension];
        
        writeImage(max(filteredStack, [], 3), xyFilteredProjectedName);
        writeImage(squeeze(max(filteredStack, [], 2)), xzFilteredProjectedName);
        writeImage(squeeze(max(filteredStack, [], 1)), yzFilteredProjectedName);
        
        if postMedian(1)
            for z = 1:size(filteredStack, 3)
                filteredStack(:, :, z) = medfilt2(filteredStack(:, :, z), [postMedian(2) postMedian(2)]);
            end;
            
            if saveStacks
                filteredStackName = [outputPath filesep '' fileNameHeader '.fusedStack.filtered_' num2str(r) '_median' outputExtension];
                writeImage(filteredStack, filteredStackName);
            end;
            
            if subProject
                firstHalfFilteredProjected  = max(filteredStack(:, :, 1:center), [], 3);
                secondHalfFilteredProjected = max(filteredStack(:, :, (center + 1):end), [], 3);
                
                firstHalfFilteredProjectedName  = [outputPath filesep '' fileNameHeader '.fusedStack_xyFProjection.filtered_' num2str(r) '_median' outputExtension];
                secondHalfFilteredProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_xyBProjection.filtered_' num2str(r) '_median' outputExtension];
                
                writeImage(firstHalfFilteredProjected, firstHalfFilteredProjectedName);
                writeImage(secondHalfFilteredProjected, secondHalfFilteredProjectedName);
            end;
            
            xyFilteredProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_xyProjection.filtered_' num2str(r) '_median' outputExtension];
            xzFilteredProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_xzProjection.filtered_' num2str(r) '_median' outputExtension];
            yzFilteredProjectedName = [outputPath filesep '' fileNameHeader '.fusedStack_yzProjection.filtered_' num2str(r) '_median' outputExtension];
            
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