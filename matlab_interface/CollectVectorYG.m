function [V_out, loop_num] = CollectVectorYG(var, data, filenames)
% [result, loop_num] = CollectVectorYG(var, data)
%
% Collect data into a vector. The function automatically loops through any
% files with an extension "*RYG.mat" in the current working directory.
%
% var : the name of the variable that should be loaded from *RYG.mat files
% data: an expression to be evaluated to generated the data to be collected
%
% For example:
%   >> var = 'Analysis'
%   >> data = 'mean(Analysis.rate{1})'
%   >> [result, loop_num] = CollectCellYG(var, data)
%
% [result, loop_num] = CollectCellYG(var, data, filenames)
%


show_progress = 0;


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

V_out = [];
loop_num = [];
fprintf('Collecting data %s from %d files: \n', data, num_files);
for i = 1:num_files
    if  show_progress == 1
        fprintf('\t Loading data %s from file %s...', data, files{i});
    end
    
    
    if strcmp(var,'ALL') %v load all
        load(files{i});
    else
        load(files{i}, var)
    end
    try
        load(files{i},  'ExplVar')
    catch
    end
    
    if  show_progress == 1
        fprintf('done.\n');
    end
    expr = sprintf('data_tmp = %s;', data);
    try        
        eval(expr);
    catch
        % warning('Cannot load %s', data)
        data_tmp = [];
    end
    
    if isempty(data_tmp)
        warning('empty data')
        files{i}
    end
    
    data_tmp = data_tmp(:)'; % row vector
    V_out = [V_out, data_tmp ]; %#ok<AGROW>
    
    if exist('ExplVar','var')
        loop_num = [loop_num, ones(1,length(data_tmp))*ExplVar.loop_num]; %#ok<AGROW>
    else
        loop_num = [];
    end
    
    clear data_tmp; % clear it! Otherwise it could be misused by the consecutive loops.
    if ~strcmp(var,'ALL') %v load all
        eval(['clear ' var]); % clear it! Otherwise it could be misused by the consecutive loops.
    end
end

fprintf('\n');
end


% function [V, loop_num] = CollectVectorYG(var)
%
% % Prepare files
% dir_strut = dir('*RYG.mat');
% num_files = length(dir_strut);
% files = cell(1,num_files);
% for id_out = 1:num_files
%     files{id_out} = dir_strut(id_out).name;
% end
%
% V = [];
% loop_num = [];
% fprintf('Collecting data %s from %d files: \n', var, num_files);
% for i = 1:num_files
%     % fprintf('Loading RYG.mat file %s...\n', files{i});
%     load(files{i});
%     if exist('R_temp','var'); % see warning in "else"
%         % disp('Loading done.');
%         eval(sprintf('data_tmp = R_temp.%s;',var));
%         data_tmp = data_tmp(:)'; % row vector
%         V = [V, data_tmp ];
%         loop_num = [loop_num, ones(1,length(data_tmp))*R_temp.ExplVar.loop_num];
%         clear R_temp; % clear it! Otherwise it could be misused by the consecutive loops.
%     else
%         warning('R_temp not found in %s! This could be due to its size being larger than 2GB.', files{i});
%     end
%     fprintf('%d,',i);
%
% end
%
% fprintf('\n');
% end