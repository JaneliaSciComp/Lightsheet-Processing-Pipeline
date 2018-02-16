function [ structure ] = convert_limits_to_values( structure, fields_to_convert )
for i=1:numel(fields_to_convert)
    structure.(fields_to_convert{i}) =  structure.(fields_to_convert{i}).start : structure.(fields_to_convert{i}).every : structure.(fields_to_convert{i}).end;
end

end

