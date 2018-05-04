% function [ output_args ] = checkInputs( input )
% fields = fieldnames(input);
% for i=1:numel(fields)
%     currentField = fields{i};
%     currentFieldValue = input.(currentField);
%     if ismember(currentField, {'inputFolder', 'outputLabel', 'inputString',...
%                                'outputString', 'outputID', 'configRoot',...
%                                'sourceString', 'inputID', 'lookUpTable',...
%                                'inputRoot','outputRoot','filePattern','inputDir', 'outputDir',...
%                                'header','footer','stackLabel'}
%         
%     elseif ismember(currentField, {})
%                    case 'thresholds'
%              validateattributes(currentFieldValue, {'double'},{'>=', 0, '<',4, 'numel', '<=', numel(input.channels)});
%     else
%         mustBeInteger(currentFieldValue)
%         switch currentField
%         
%             case 'specimen'
%                 validateattributes(currentFieldValue, {'double'}, {'>=',0,'<',100})
%             case {'timepoints'}
%             
%             
%     end
%      switch currentField
%         case 'specimen'
%             %validateattributes(currentFieldValue, {'double'}, {'>=',0,'<',100})
%             mustBeInteger(currentFieldValue);
%        % case {'timepoints'}
%        %     mustBeInteger...
%         case {'cameras', 'channels'}
%             validateattributes(currentFieldValue, {'double'},{'<',4, 'numel', '<=',4})
%         case 'dimensions'
%           %  mustBeInteger(currentFieldValue);
%         case {'startsLeft', 'startsTop', 'widths', 'heights', 'startsFront', 'depths'}
%             if ~isempty(currentFieldValue)            
%                 validateattributes(currentFieldValue, {'double'}, {'numel','=', numel(input.cameras)})
%             end
%          case 'segmentFlag'
%              validateattributes(currentFieldValue, {'double'},{'>=', 0, '<',4});
%         
%          %case references dependents
% 
%     end
%     end
% 
% end
% function checkIfValidFile
% if ~isempty(regexp(filename, '[/\*:?"<>|]', 'once'))
%         %Do something
% end