function PostPostProcessYG(stdin)
% This function is not very useful

 % Prepare files
 if nargin == 0
%      wd = cd; % store current directory
%      cd /import/yossarian1/yifan/Project1/
%      addpath(genpath(cd));
%      cd(wd); % return
     dir_strut = dir('*RYG.mat');
     num_files = length(dir_strut);
     files = cell(1,num_files);
     for id_out = 1:num_files
         files{id_out} = dir_strut(id_out).name;
     end
 else
     % stdin, i.e., file pathes and names separated by space
     files = textscan(stdin,'%s'); % cell array of file path+names
     num_files = length(files);
     for i = 1:num_files
         files{i} = cell2mat(files{i});
     end
 end
 
 % Post-postprocess data
 save_fig = -1;
 for i = 1:num_files
     fprintf('Loading RYG.mat file %s...\n', files{i});
     R_temp = load(files{i},'reduced','N');
     
     disp('Loading done.');
     %%%%%%% do something here
     % R = ClusterYG({R_temp},save_fig);
     % R_temp = cluster_sorted_rate(R_temp);
     % R_temp = R{1};
     R_temp = get_stPR(R_temp);
     save(files{i},'-struct', 'R_temp', '-v7.3','-append'); % -v7.3 for >2GB
     %%%%%%% visualise results
     % RasterYG(R, save_fig);
     % ClusterYG(R,save_fig);
     % HistogramsYG(R,save_fig);
     
 end

% var = 'Analysis.Hz_overall';
% [Hz, loop_num] = CollectVectorYG(var);
% var = 'ExplVar.EE_factor';
%[EE, dummy] = CollectVectorYG(var);
% var = 'ExplVar.II_factor';
% [II, dummy] = CollectVectorYG(var);
 
% save('balance_Hz','Hz','loop_num','EE','II');


% var = 'up_down_analysis.up_C_label';
% [up_C_label, loop_num_up_C_label] = CollectVectorYG(var);
% var = 'up_down_analysis.up_duration';
% [up_duration, dummy] = CollectVectorYG(var);
% var = 'up_down_analysis.down_duration';
% [down_duration, dummy] = CollectVectorYG(var);
% var = 'up_down_analysis.up_overlapping';
% [up_overlapping, dummy] = CollectVectorYG(var);
% var = 'Analysis.Hz_overall';
% [Hz, loop_num_Hz] = CollectVectorYG(var);
% save('phase_diagram_data_2','up_C_label','up_duration','down_duration',...
% 	'up_overlapping','Hz','loop_num_up_C_label','loop_num_Hz');

% 
% var = 'cluster.low_du{1}';
% [low_du, loop_num] = CollectVectorYG(var);
%  
% save('low_du','low_du','loop_num');

end
