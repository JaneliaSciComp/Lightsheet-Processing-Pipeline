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
str = sprintf('mcc -m -R -nodesktop -v clusterPT_fn.m %s;', astr);
eval(str);
