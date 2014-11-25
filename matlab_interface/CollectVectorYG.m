function [V, loop_num] = CollectVectorYG(var, data)
% example:
% var = 'cluster'
% data = 'cluster.high_du{3}'


% Prepare files
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end

V = [];
loop_num = [];
fprintf('Collecting data %s from %d files: \n', data, num_files);
for i = 1:num_files
    fprintf('\t Loading data %s from file %s...', data, files{i});
    load(files{i}, var, 'ExplVar');

    fprintf('done.\n');
    
        eval(sprintf('data_tmp = %s;', data));
        data_tmp = data_tmp(:)'; % row vector
        V = [V, data_tmp ];
        loop_num = [loop_num, ones(1,length(data_tmp))*ExplVar.loop_num];
        clear data_tmp; % clear it! Otherwise it could be misused by the consecutive loops.
        
    
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