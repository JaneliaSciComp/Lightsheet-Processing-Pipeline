inputName = ['R:' filesep 'SV1' filesep 'KM_16-03-28' filesep 'Mmu_E1_H2BmCherry_20160328_170713.corrected' filesep 'Results' filesep 'TimeFused.Corrected' filesep 'Mmu_H2BmCherry.TM000050_timeFused_blending' filesep 'SPM00_TM000050_CM00_CM01_CHN00.fusedStack.corrected.klb'];

% copyfile(inputName, [inputName(1:(end-3)) 'old.klb']);
% copyfile([inputName(1:(end-14)) '_xyProjection.corrected.klb'], [inputName(1:(end-14)) '_xyProjection.corrected.old.klb']);
% copyfile([inputName(1:(end-14)) '_xzProjection.corrected.klb'], [inputName(1:(end-14)) '_xzProjection.corrected.old.klb']);
% copyfile([inputName(1:(end-14)) '_yzProjection.corrected.klb'], [inputName(1:(end-14)) '_yzProjection.corrected.old.klb']);

stack = readImage(inputName);
stack(:,1:(end-3),247) = stack(:,4:end,247);
stack(16:end,:,1:246) = stack(1:(end-15),:,1:246);

writeImage(stack, inputName);
writeImage(max(stack,[],3), [inputName(1:(end-14)) '_xyProjection.corrected.klb']);
writeImage(squeeze(max(stack,[],2)), [inputName(1:(end-14)) '_xzProjection.corrected.klb']);
writeImage(squeeze(max(stack,[],1)), [inputName(1:(end-14)) '_yzProjection.corrected.klb']);