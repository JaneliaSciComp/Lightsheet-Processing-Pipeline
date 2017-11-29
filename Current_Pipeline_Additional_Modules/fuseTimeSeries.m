function [ output_args ] = fuseTimeSeries( directory,  file_suffix)
all_directories = dir(directory);
all_directories = all_directories(3:end);
all_stacks = cell(numel(all_directories),1);
dimensions = zeros(numel(all_directories),3);

parfor current_directory_index = 1:numel(all_directories)
    directory_path = [all_directories(current_directory_index).folder '/' all_directories(current_directory_index).name '/' ];
    file_name = dir([directory_path '*' file_suffix]);
    current_stack = readKLBstack([directory_path file_name.name]);
    all_stacks{current_directory_index} = current_stack;
    dimensions(current_directory_index,:) = size(current_stack);
end

fused_stack = zeros([max(dimensions), 1, numel(all_directories)],'uint8');
max_z = max(dimensions(:,3));
for current_time_point = 1:numel(all_directories)
  %  current_stack = zeros(max(dimensions),'uint8');
    current_z_size = dimensions(current_time_point,3);
  %  start_z = (max_z-current_z_size)/2 + 1;
    %current_stack(:,:, start_z : start_z + current_z_size-1) = uint8(all_stacks{current_time_point});
    %    current_stack(:,:, 1 : current_z_size) = uint8(all_stacks{current_time_point});
    fused_stack(:,:, 1:current_z_size, :, current_time_point) = uint8(all_stacks{current_time_point});
end
end

