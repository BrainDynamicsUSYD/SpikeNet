function [ R ] = CollectRYG(varargin)
%This function collects R_temp data in the default path data/
% input should be dir(.)
%   Detailed explanation goes here
%   This function is not very useful due to the large memeory requirement

disp('Collecting RYG files...');
tic;

% Prepare filenames, "name" should be cell array of filename strings (*.ygout)
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

% Collecting data
R = cell(1,num_files);
for i = 1:num_files
    load(strcat('data/',name{i}));
    if any( strcmp(fieldnames(R_temp.ExplVar), 'loop_num') )
        loop_num = R_temp.ExplVar.loop_num;
    else
        loop_num = i;
    end
    R{loop_num} = R_temp;
end

toc;
disp('Collecting done.');


% % from cell array of struct to struct array, only making sense when the
% % fields are identical!
% Result_cell = cell2mat(Result_cell); 

end

