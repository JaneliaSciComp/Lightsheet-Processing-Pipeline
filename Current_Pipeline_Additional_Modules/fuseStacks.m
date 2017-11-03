timepoints       = 0:278;
dataTypes        = {'LargePSF_iter50'; 'SmallPSF_iter20'};
sigmas           = [20 8; 5 2];
kernelFactor     = 3;
splitting        = 5;
writeProjections = 1;

if exist('ProjectionsFusion', 'dir') ~= 7
    mkdir('ProjectionsFusion');
end;

if exist('ProjectionsFusion' filesep 'Frames', 'dir') ~= 7
    mkdir('ProjectionsFusion' filesep 'Frames');
end;

for t = 1:numel(timepoints)
    timepoint = timepoints(t);
    
    for d = 1:numel(dataTypes)
        sigma1 = sigmas(d, 1);
        sigma2 = sigmas(d, 2);
        
        outputStackName = ['TM' num2str(timepoint, '%.6d') filesep 'SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00' ...
            '.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_' dataTypes{d} '_lambdaTV000000.fusionSigma_' num2str(sigma1) '_' num2str(sigma2) '.uint16.klb'];
        
        lastProjectionName = ['ProjectionsFusion' filesep 'Frames' filesep 'Fusion_' dataTypes{d} '.fusionSigma_' num2str(sigma1) '_' num2str(sigma2) '.TM' num2str(timepoint, '%.3d') '_yz.klb'];
        
        if exist(outputStackName, 'file') ~= 2 || (writeProjections && exist(lastProjectionName, 'file') ~= 2)
            disp(' ');
            disp(['Reading data type ' dataTypes{d} ' at time point ' num2str(timepoint)]);
            
            stacks = cell(4, 1);
            
            for s = 1:4
                switch s
                    case 1
                        specimenString = 'SPM00';
                        cameraString = 'CM00';
                    case 2
                        specimenString = 'SPM00';
                        cameraString = 'CM01';
                    case 3
                        specimenString = 'SPM01';
                        cameraString = 'CM00';
                    case 4
                        specimenString = 'SPM01';
                        cameraString = 'CM01';
                end;
                
                stacks{s, 1} = readImage(['TM' num2str(timepoint, '%.6d') filesep '' specimenString '_TM' num2str(timepoint, '%.6d') '_' cameraString '_CHN00' ...
                    '.affine.trsf.cropped.klb_dec_LR_multiGPU_SVD_' dataTypes{d} '_lambdaTV000000.uint16.klb']);
            end;
            
            kernelSize1  = [1 1 1] .* (sigma1*kernelFactor);
            kernelSigma1 = [1 1 1] .* sigma1;
            
            kernelSize2  = [1 1 1] .* (sigma2*kernelFactor);
            kernelSigma2 = [1 1 1] .* sigma2;
            
            weightingFactors = cell(4, 1);
            
            for i = 1:4
                disp(['Creating weighting factors for stack ' num2str(i) ' using sigma1 = ' num2str(sigma1) ' and sigma2 = ' num2str(sigma2)]);
                
                if splitting > 1
                    splittingMargin = 2 * sigma1 * kernelFactor;
                    processedStack = zeros(size(stacks{i}, 1), size(stacks{i}, 2), size(stacks{i}, 3));
                    for s = 1:splitting
                        xStart = max(1, round((s - 1) * size(stacks{i}, 1) / splitting + 1 - splittingMargin));
                        xStop = min(size(stacks{i}, 1), round(s * size(stacks{i}, 1) / splitting + splittingMargin));
                        convolvedSlab = imgaussianAnisotropy(double(stacks{i}(xStart:xStop, :, :)), kernelSigma1, kernelSize1);
                        if s == 1
                            processedStack(1:(xStop - splittingMargin), :, :) = convolvedSlab(1:(end - splittingMargin), :, :);
                        elseif s == splitting
                            processedStack((xStart + splittingMargin):end, :, :) = convolvedSlab((1 + splittingMargin):end, :, :);
                        else % s > 1 && s < splitting
                            processedStack((xStart + splittingMargin):(xStop - splittingMargin), :, :) = convolvedSlab((1 + splittingMargin):(end - splittingMargin), :, :);
                        end;
                        clear convolvedSlab;
                    end;
                    processedStack = (double(stacks{i}) - processedStack) .^ 2;
                    
                    splittingMargin = 2 * sigma2 * kernelFactor;
                    weightingFactors{i} = zeros(size(stacks{i}, 1), size(stacks{i}, 2), size(stacks{i}, 3));
                    for s = 1:splitting
                        xStart = max(1, round((s - 1) * size(stacks{i}, 1) / splitting + 1 - splittingMargin));
                        xStop = min(size(stacks{i}, 1), round(s * size(stacks{i}, 1) / splitting + splittingMargin));
                        convolvedSlab = imgaussianAnisotropy(processedStack(xStart:xStop, :, :), kernelSigma2, kernelSize2);
                        if s == 1
                            weightingFactors{i}(1:(xStop - splittingMargin), :, :) = convolvedSlab(1:(end - splittingMargin), :, :);
                        elseif s == splitting
                            weightingFactors{i}((xStart + splittingMargin):end, :, :) = convolvedSlab((1 + splittingMargin):end, :, :);
                        else % s > 1 && s < splitting
                            weightingFactors{i}((xStart + splittingMargin):(xStop - splittingMargin), :, :) = convolvedSlab((1 + splittingMargin):(end - splittingMargin), :, :);
                        end;
                        clear convolvedSlab;
                    end;
                    clear processedStack;
                else
                    weightingFactors{i} = (double(stacks{i}) - imgaussianAnisotropy(double(stacks{i}), kernelSigma1, kernelSize1)) .^ 2;
                    weightingFactors{i} = imgaussianAnisotropy(weightingFactors{i}, kernelSigma2, kernelSize2);
                end;
            end;
            
            disp(['Fusing stacks for sigma1 = ' num2str(sigma1) ' and sigma2 = ' num2str(sigma2)]);
            
            for i = 1:4
                if i == 1
                    fusedStack = weightingFactors{i} .* double(stacks{i});
                    factorSum = double(weightingFactors{i});
                    weightingFactors{i} = [];
                    stacks{i} = [];
                else
                    fusedStack = fusedStack + weightingFactors{i} .* double(stacks{i});
                    factorSum = factorSum + double(weightingFactors{i});
                    weightingFactors{i} = [];
                    stacks{i} = [];
                end;
            end;
            
            fusedStack = uint16(fusedStack ./ factorSum);
            
            clear factorSum;
            
            writeImage(fusedStack, outputStackName);
            
            if writeProjections
                writeImage(max(fusedStack, [], 3), ['ProjectionsFusion' filesep 'Frames' filesep 'Fusion_' dataTypes{d} '.fusionSigma_' num2str(sigma1) '_' num2str(sigma2) '.TM' num2str(timepoint, '%.3d') '_xy.klb']);
                writeImage(squeeze(max(fusedStack, [], 2)), ['ProjectionsFusion' filesep 'Frames' filesep 'Fusion_' dataTypes{d} '.fusionSigma_' num2str(sigma1) '_' num2str(sigma2) '.TM' num2str(timepoint, '%.3d') '_xz.klb']);
                writeImage(squeeze(max(fusedStack, [], 1)), ['ProjectionsFusion' filesep 'Frames' filesep 'Fusion_' dataTypes{d} '.fusionSigma_' num2str(sigma1) '_' num2str(sigma2) '.TM' num2str(timepoint, '%.3d') '_yz.klb']);
            end;
            
            clear fusedStack;
        end;
    end;
end;

disp(' ');