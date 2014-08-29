function [V, loop_num] = CollectVectorYG(var)

% Prepare files
dir_strut = dir('*RYG.mat');
num_files = length(dir_strut);
files = cell(1,num_files);
for id_out = 1:num_files
    files{id_out} = dir_strut(id_out).name;
end

% Post-postprocess data
save_fig = 1;
V = [];
loop_num = [];
for i = 1:num_files
    % fprintf('Loading RYG.mat file %s...\n', files{i});
    load(files{i}); 
    % disp('Loading done.');
    %%%%%%% collect data here
    eval(sprintf('data_tmp = R_temp.%s;',var));
    data_tmp = data_tmp(:)'; % row vector
    V = [V, data_tmp ];
    loop_num = [loop_num, ones(1,length(data_tmp))*R_temp.ExplVar.loop_num];
    %%%%%%% visualise results
    % RasterYG(R, save_fig);
    % ClusterYG(R,save_fig);
    % HistogramsYG(R,save_fig);
    
end

end