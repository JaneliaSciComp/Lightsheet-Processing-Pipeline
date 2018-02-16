function compile_function(function_name)
%% compile (you must customize this script to your system/environment


cd /groups/lightsheet/lightsheet/home/ackermand/Lightsheet-Processing-Pipeline/Compiled_Functions/


directories = {...
    '/groups/lightsheet/lightsheet/home/ackermand/Lightsheet-Processing-Pipeline/Current_Pipeline/',...
    '/groups/lightsheet/lightsheet/home/ackermand/Lightsheet-Processing-Pipeline/Current_Pipeline_Additional_Modules/',...
    '/groups/lightsheet/lightsheet/home/ackermand/Lightsheet-Processing-Pipeline/Pipeline_Functions/',...
    '/groups/lightsheet/lightsheet/home/ackermand/Lightsheet-Processing-Pipeline/external/jsonlab/',...
    };

astr = [];
for current_dir=1:numel(directories)
    fn = dir([directories{current_dir} '*.m']);
    for ix = 1:numel(fn)
        astr = [astr sprintf(' -a %s%s',directories{current_dir},fn(ix).name)];
    end
end
directory = '/groups/lightsheet/lightsheet/home/ackermand/compile_code/keller-lab-block-filetype/matlabWrapper/';
fn = dir([directory '*.mexa64']);
for ix = 1:numel(fn)
    astr = [astr sprintf(' -a %s%s',directory,fn(ix).name)];
end
if isequal(function_name,'localEC_fn')
   str=sprintf('mcc -m -R -nodesktop -v -w off:MATLAB:mir_warning_maybe_uninitialized_temporary %s.m %s;', function_name, astr) %Had to add this because it produces a warning otherwise, and within that warning is the word "error" which I think might be a problem in how JACS determines if an error occurs
else
    str = sprintf('mcc -m -R -nodesktop -v %s.m %s;', function_name, astr);
end
eval(str);
end