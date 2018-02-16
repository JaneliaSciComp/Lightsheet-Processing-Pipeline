function distribute_parallelized( function_name, filename, processors_per_node, timepoints_per_node, waitfor_job_name )
if nargin==2
    timepoints_per_node = processors_per_machine;
end
if nargin==4
   waitfor_job_name = []; 
end
input_parameters = loadjson(fileread(filename));
if ~isfield(input_parameters, 'verbose')
    input_parameters.verbose = true;
end
if input_parameters.verbose == true
    disp(['Using input file: ' filename]);
    disp(input_parameters)
end
if ischar(processors_per_node), processors_per_node = str2double(processors_per_node); end
if ischar(timepoints_per_node), timepoints_per_node = str2double(timepoints_per_node); end
number_of_nodes = ceil(numel(input_parameters.timepoints)/timepoints_per_node);
timepoints = input_parameters.timepoints;
fprintf('Running %d timepoints on %d nodes\n',numel(input_parameters.timepoints), number_of_nodes);
%for job_number=1:number_of_nodes
  %  input_parameters.timepoints = timepoints(~isnan(timepoints(:,job_number)),job_number);
%    system(sprintf('bsub -P lightsheet -n %d ./run_compiled_matlab.sh ./run_%s_fn.sh /misc/local/matlab-2017a/ %s %d %d &', processors_per_node, function_name, filename, timepoints_per_node, job_number));
%end
if contains(function_name, 'local')
     if isempty(waitfor_job_name) || waitfor_job_name=="none"
      %  system(sprintf('bsub -J "%s" -P lightsheet -n %d /groups/lightsheet/lightsheet/home/ackermand/Lightsheet-Processing-Pipeline/Compiled_Functions/run_compiled_matlab.sh /groups/lightsheet/lightsheet/home/ackermand/Lightsheet-Processing-Pipeline/Compiled_Functions/run_%s_fn.sh /misc/local/matlab-2017a/ %s', function_name, processors_per_node, function_name));
    else
     %   system(sprintf('bsub -w "done(%s)" -ti -J "%s" -P lightsheet -n %d /groups/lightsheet/lightsheet/home/ackermand/Lightsheet-Processing-Pipeline/Compiled_Functions/run_compiled_matlab.sh /groups/lightsheet/lightsheet/home/ackermand/Lightsheet-Processing-Pipeline/Compiled_Functions/run_%s_fn.sh /misc/local/matlab-2017a/ %s', waitfor_job_name, function_name, processors_per_node, function_name, filename));
    end
else
    if isempty(waitfor_job_name) || waitfor_job_name=="none"
     %   system(sprintf('bsub -J "%s[1-%d]" -P lightsheet -n %d /groups/lightsheet/lightsheet/home/ackermand/Lightsheet-Processing-Pipeline/Compiled_Functions/run_compiled_matlab.sh /groups/lightsheet/lightsheet/home/ackermand/Lightsheet-Processing-Pipeline/Compiled_Functions/run_%s_fn.sh /misc/local/matlab-2017a/ %s %d \\$LSB_JOBINDEX', function_name, number_of_nodes, processors_per_node, function_name, filename, timepoints_per_node));
    else
     %   system(sprintf('bsub -w "done(%s)" -ti -J "%s[1-%d]" -P lightsheet -n %d /groups/lightsheet/lightsheet/home/ackermand/Lightsheet-Processing-Pipeline/Compiled_Functions/run_compiled_matlab.sh /groups/lightsheet/lightsheet/home/ackermand/Lightsheet-Processing-Pipeline/Compiled_Functions/run_%s_fn.sh /misc/local/matlab-2017a/ %s %d \\$LSB_JOBINDEX', waitfor_job_name, function_name, number_of_nodes, processors_per_node, function_name, filename, timepoints_per_node));
    end
end

