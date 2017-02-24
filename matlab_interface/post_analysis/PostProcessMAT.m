function PostProcessMAT(stdin)
% stdin is paths of input files in unix-style, i.e., separated by spaces
% If given no argument, it searches for matches under CURRENT directory


% Prepare files
if nargin == 0
    dir_strut = dir('*RYG.mat');
    num_files = length(dir_strut);
    files = cell(1,num_files);
    for i = 1:num_files
        files{i} = dir_strut(i).name;
    end
else
    % stdin, i.e., file pathes and names separated by space
    files = textscan(stdin,'%s'); % cell array of file path+names
    num_files = length(files);
    for i = 1:num_files
        files{i} = cell2mat(files{i});
    end
end

% save figures
save_fig = 1; % -1 for no figure, 0 for displaying figure, 1 for saving figure
% Start processing
for i = 1:num_files
    
    % start form .mat files
    fprintf('Loading RYG.mat file %s...', files{i});
    R_temp = load(files{i}); % for performance, only load the necessary stuff
    disp('done.');
    %%%%%%% do something here


    
    
%     R_temp = avalanche_detect(R_temp);
%     avalanche = R_temp.avalanche;
%     save(files{i},'avalanche', '-append');
%         
%     R_temp = get_stPR(R_temp);
%     stPR = R_temp.stPR;
%     save(files{i},'stPR', '-append');
%     

    

%      R_temp = get_SWR(R_temp);
%      LFP = R_temp.LFP;
%      save(files{i},'LFP', '-append');
%      
%      [R_temp] = get_CC_pop(R_temp, 1);
%      Analysis = R_temp.Analysis;
%      save(files{i},'Analysis', '-append');
%      

%      [R_temp] = get_lagged_cov(R_temp);
%      Analysis = R_temp.Analysis;
%      save(files{i},'Analysis', '-append');

      get_LFP_continous(R_temp);


          R_temp = get_grid_firing_centre(R_temp,'mode','bayesian');
          grid = R_temp.grid; %#ok<NASGU>
          save(files{i},'grid', '-append');
     
         R_temp = get_grid_SWR_consistency(R_temp);
         grid_SWR = R_temp.grid_SWR; %#ok<NASGU>
         save(files{i},'grid_SWR', '-append');
     
     R_temp = get_SWR_spike_phase_lock(R_temp);
     SWR_spike_phase_lock = R_temp.SWR_spike_phase_lock; %#ok<NASGU>
     save(files{i},'SWR_spike_phase_lock', '-append');
     
     % plot_SWR(R_temp, save_fig);
     %
    
    
    %     R_temp = get_CC_pop(R_temp);
    %     R_temp = get_EI_current_crosscorr(R_temp);
    %     Balance = R_temp.Balance;
    %     Analysis = R_temp.Analysis;
    %     save(files{i},'Balance', 'Analysis', '-append');
    

    
    %save(files{i},'-struct', 'R_temp', '-v7.3'); % -v7.3 for >2GB
    
    %     R_temp = get_neuron_sample_stats(R_temp);
    %     neuron_sample_stats = R_temp.neuron_sample_stats;
    %     save(files{i},'neuron_sample_stats', '-append');
    
    
    
    %%%%%%% do something above
    disp('Data processed and saved.');
    
end


end
