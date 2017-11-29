%% compile (you must customize this script to your system/environment


cd /groups/lightsheet/lightsheet/home/ackermand/Lightsheet-Processing-Pipeline/Compiled_Functions/


directories = {...
    '/groups/lightsheet/lightsheet/home/ackermand/Lightsheet-Processing-Pipeline/external/jsonlab/',...
    };

astr = [];
for current_dir=1:numel(directories)
    fn = dir([directories{current_dir} '*.m']);
    for ix = 1:numel(fn)
        astr = [astr sprintf(' -a %s%s',directories{current_dir},fn(ix).name)];
    end
end
str = sprintf('mcc -m -R -nodesktop -v distribute_parallelized.m %s;', astr);
eval(str);
