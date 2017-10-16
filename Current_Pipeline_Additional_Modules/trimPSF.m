threshold = 1e-10;

candidates = dir;
for i = numel(candidates):-1:1
    if ~isdir(candidates(i).name) || strcmp(candidates(i).name, '.') || strcmp(candidates(i).name, '..') || strcmp(candidates(i).name, 'XMLs')
        candidates(i) = [];
    else
        timepoint = str2num(candidates(i).name(3:end));
        missingFlag = ...
            exist([candidates(i).name '\PSF\PSF_YW_2.2_5_SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\PSF\PSF_YW_2.2_5_SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\PSF\PSF_YW_2.2_5_SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\PSF\PSF_YW_2.2_5_SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\PSF\PSF_ZZ_2.5_8_SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\PSF\PSF_ZZ_2.5_8_SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\PSF\PSF_ZZ_2.5_8_SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.klb'], 'file') ~= 2 || ...
            exist([candidates(i).name '\PSF\PSF_ZZ_2.5_8_SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.klb'], 'file') ~= 2;
        if missingFlag
            candidates(i) = [];
        else
            alldoneFlag = ...
                exist([candidates(i).name '\PSF\PSF_YW_2.2_5_SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.trimmed.klb'], 'file') == 2 && ...
                exist([candidates(i).name '\PSF\PSF_YW_2.2_5_SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.trimmed.klb'], 'file') == 2 && ...
                exist([candidates(i).name '\PSF\PSF_YW_2.2_5_SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.trimmed.klb'], 'file') == 2 && ...
                exist([candidates(i).name '\PSF\PSF_YW_2.2_5_SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.trimmed.klb'], 'file') == 2 && ...
                exist([candidates(i).name '\PSF\PSF_ZZ_2.5_8_SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.trimmed.klb'], 'file') == 2 && ...
                exist([candidates(i).name '\PSF\PSF_ZZ_2.5_8_SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.trimmed.klb'], 'file') == 2 && ...
                exist([candidates(i).name '\PSF\PSF_ZZ_2.5_8_SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.trimmed.klb'], 'file') == 2 && ...
                exist([candidates(i).name '\PSF\PSF_ZZ_2.5_8_SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.trimmed.klb'], 'file') == 2;
            if alldoneFlag
                candidates(i) = [];
            end;
        end;
    end;
end;

for i = 1:numel(candidates)
    folder = [candidates(i).name '\PSF\'];
    timepoint = str2num(candidates(i).name(3:end));
    
    disp(['trimming PSF at timepoint ' num2str(timepoint)]);
    
    for f = 1:8
        switch f
            case 1
                psfName = [folder 'PSF_YW_2.2_5_SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.klb'];
            case 2
                psfName = [folder 'PSF_YW_2.2_5_SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.klb'];
            case 3
                psfName = [folder 'PSF_YW_2.2_5_SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.klb'];
            case 4
                psfName = [folder 'PSF_YW_2.2_5_SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.klb'];
            case 5
                psfName = [folder 'PSF_ZZ_2.5_8_SPM00_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.klb'];
            case 6
                psfName = [folder 'PSF_ZZ_2.5_8_SPM00_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.klb'];
            case 7
                psfName = [folder 'PSF_ZZ_2.5_8_SPM01_TM' num2str(timepoint, '%.6d') '_CM00_CHN00.klb'];
            case 8
                psfName = [folder 'PSF_ZZ_2.5_8_SPM01_TM' num2str(timepoint, '%.6d') '_CM01_CHN00.klb'];
        end;
        
        PSF = readImage(psfName);
        PSF = PSF ./ sum(PSF(:));
        
        BW = PSF > (max(PSF(:)) * threshold);
        
        erase = zeros(3, 1);
        
        for ii = 1:ndims(BW)
            qq = BW;
            for jj = ndims(BW):-1:1
                if ii == jj
                    continue;
                end;
                qq = sum(qq, jj);
            end;
            qq = squeeze(qq);
            
            p1 = find(qq > 0, 1, 'first') - 2;
            p1 = max(1, p1);
            p2 = find(qq > 0, 1, 'last') + 2;
            p2 = min(length(qq), p2);
            
            erase(ii) = min(p1, length(qq(p2:end)));
        end;
        
        PSF = PSF(erase(1):(end - erase(1) + 1), erase(2):(end - erase(2) + 1), erase(3):(end - erase(3) + 1));
        
        PSF = PSF ./ sum(PSF(:));
        
        writeImage(single(PSF), [psfName(1:(end - 3)) 'trimmed' psfName((end - 3):end)]);
    end;
end;