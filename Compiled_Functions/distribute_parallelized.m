function distribute_parallelized( function_name, filename, processors_per_node, timepoints_per_node )
if nargin==2
    timepoints_per_node = processors_per_machine;
end
input_parameters = loadjson(fileread(filename));
if input_parameters.verbose
    disp(['Using input file: ' filename]);
    disp(input_parameters)
end
if ischar(processors_per_node), processors_per_node = str2double(processors_per_node); end
if ischar(timepoints_per_node), timepoints_per_node = str2double(timepoints_per_node); end
number_of_nodes = ceil(numel(input_parameters.timepoints)/timepoints_per_node);
timepoints = input_parameters.timepoints;
fprintf('Running %d timepoints on %d nodes\n',numel(input_parameters.timepoints), number_of_nodes);
for job_number=1:number_of_nodes
    input_parameters.timepoints = timepoints(~isnan(timepoints(:,job_number)),job_number);
    system(sprintf('bsub -P lightsheet -n %d ./run_compiled_matlab.sh ./run_%s_fn.sh /misc/local/matlab-2017a/ %s %d %d &', processors_per_node, function_name, filename, timepoints_per_node, job_number));
end

