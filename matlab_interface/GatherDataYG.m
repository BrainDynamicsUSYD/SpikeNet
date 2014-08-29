function GatherDataYG(varargin)

%%%%%%%%%%% Prepare filenames
% "name" should be cell array of filename strings (*.ygout)
if ~isempty(varargin)
    dir_strut = varargin{1};
else % read all
    dir_strut = dir('data/*RYG.mat');
end
num_files = length(dir_strut);
name = cell(1,num_files);
for id_out = 1:num_files
    name{id_out} = dir_strut(id_out).name;
end


%%%%%%%%%%%% gather data-of-interest
ExplVar = cell(1,num_files);
up_down = cell(1,num_files);
up_duration = [];
down_duration = [];
for i = 1:num_files
    % load raw data
    load(strcat('data/',name{i}));

    % gather data-of-interest
    up_down{i} = R_temp.up_down_analysis;
    up_duration = [up_duration up_down{i}.up_duration];
    down_duration = [down_duration up_down{i}.down_duration];
    
    ExplVar{i} = R_temp.ExplVar;
    disp(i/num_files)
end
up_down = cell2mat(up_down);
ExplVar = cell2mat(ExplVar);




%%%%%%%%%%%%% save data-of-interest
save('/import/yossarian1/yifan/Project1/up_down_analysis.mat', 'up_down','up_duration','down_duration', 'ExplVar');

end

