function [V_tmp, loop_num] = CollectCellYG(var, data, filenames)
% [result, loop_num] = CollectCellYG(var, data)
%
% Collect data into a cell. The function automatically loops through any
% files with an extension "*RYG.mat" in the current working directory.
% 
% var : the name of the variable that should be loaded from *RYG.mat files
% data: an expression to be evaluated to generated the data to be collected
%
% For example:
%   >> var = 'cluster'
%   >> data = 'mean(cluster.high_du{3})'
%   >> [result, loop_num] = CollectCellYG(var, data)
% 
% [result, loop_num] = CollectCellYG(var, data, filenames)
%

% Prepare files
if nargin == 2
    dir_strut = dir('*RYG.mat');
elseif nargin == 3
    dir_strut = dir(filenames);
end

num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end

V_tmp = cell(1,num_files);
loop_num = [];
fprintf('Collecting data %s from %d files: \n', data, num_files);
for i = 1:num_files
    % fprintf('\t Loading data %s from file %s...', data, files{i});
    if strcmp(var,'ALL') %v load all
        load(files{i});
    else
        load(files{i}, var)
    end
    try 
        load(files{i},  ExplVar')
    catch
    end
    
    data_tmp = [];
    try        
        eval(sprintf('data_tmp = %s;', data));
    catch
        % warning('Cannot load %s', data)
    end
    V_tmp{i} = data_tmp;
    
    if exist('ExplVar','var')
        loop_num = [loop_num, ExplVar.loop_num]; %#ok<AGROW>
    else
        loop_num = [];
    end

    clear data_tmp; % clear it! Otherwise it could be misused by the consecutive loops.
    if ~strcmp(var,'ALL') %v load all
        eval(['clear ' var]); % clear it! Otherwise it could be misused by the consecutive loops.
    end
    %fprintf('done.\n');
end

fprintf('\n');
end

