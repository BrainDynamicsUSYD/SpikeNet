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

%     R_temp = get_CN_prob(R_temp);
%     CN_prob = R_temp.CN_prob;
%     save(files{i},'CN_prob', '-append');
    
%     R_temp = avalanche_detect(R_temp);
%     avalanche = R_temp.avalanche;
%     save(files{i},'avalanche', '-append');
        
%     R_temp = get_stPR(R_temp,'n_sample_region', 2000);
%     stPR = R_temp.stPR;
%     save(files{i},'stPR', '-append');
%     
% 
%     
% 
% %      R_temp = get_SWR(R_temp);
% %      LFP = R_temp.LFP;
% %      save(files{i},'LFP', '-append');
% %      
     [R_temp] = get_CC_pop(R_temp, 1);
     % [R_temp] = get_dist_CC(R_temp);
     Analysis = R_temp.Analysis;
     save(files{i},'Analysis', '-append');
% %      

%      [R_temp] = get_lagged_cov(R_temp);
%      Analysis = R_temp.Analysis;
%      save(files{i},'Analysis', '-append');

%       get_LFP_continous(R_temp);


%            R_temp = get_grid_firing_centre(R_temp,'mode','bayesian');
%            grid = R_temp.grid; %#ok<NASGU>
%            grid_sub = R_temp.grid_sub; %#ok<NASGU>
%            save(files{i},'grid','grid_sub', '-append');

%            R_temp = get_grid_firing_centre(R_temp,'mode','quick','win_len', 150);
%            grid_150 = R_temp.grid; %#ok<NASGU>
%            save(files{i},'grid_150', '-append');

%         hw = 31;
% % [R_temp] = get_fano_factor(R_temp, 2, hw);
% [R_temp] = get_fano_factor(R_temp, 1);
% Analysis = R_temp.Analysis;
% save(files{i},'Analysis', '-append');

% 
%[R_temp] = get_rich_club(R_temp);
%rich_club = R_temp.rich_club;
% save(files{i},'rich_club', '-append');

% R_temp = get_triplet_sequence(R_temp);
% triplet = R_temp.triplet;
% save(files{i},'triplet', '-append');

% R_temp = get_triplet_sequence(R_temp,'hw_sample', 10);
% triplet_double = R_temp.triplet;
% save(files{i},'triplet_double', '-append');


% R_temp = get_triplet_sequence_no_latency(R_temp,'n_trial', 50);
% triplet_no_latency = R_temp.triplet;
% save(files{i},'triplet_no_latency', '-append');

% R_temp = get_local_spike_corr(R_temp);
% local_cc = R_temp.local_cc;
% save(files{i},'local_cc', '-append');

% % 
% % 
% [R_temp] = get_motif(R_temp);
% motif = R_temp.motif;
% save(files{i},'motif', '-append');

% 
% R_temp = get_stru_var_decomp(R_temp);
% stru_var_decomp = R_temp.stru_var_decomp;
% save(files{i},'stru_var_decomp', '-append');


%          R_temp = get_grid_SWR_consistency(R_temp);
%          grid_SWR = R_temp.grid_SWR; %#ok<NASGU>
%          save(files{i},'grid_SWR', '-append');
% 
%          R_temp = get_grid_SWR_consistency(R_temp);
%          grid_SWR_6 = R_temp.grid_SWR; %#ok<NASGU>
%          save(files{i},'grid_SWR_6', '-append');

%      R_temp = get_SWR_spike_phase_lock(R_temp);
%      SWR_spike_phase_lock = R_temp.SWR_spike_phase_lock; %#ok<NASGU>
%      save(files{i},'SWR_spike_phase_lock', '-append');

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
    

% R_temp = get_weight_vs_CC(  R_temp );
% w_CC = R_temp.w_CC;
% save(files{i},'w_CC', '-append');
    %%%%%%% do something above
    disp('Data processed and saved.');
    
end


end
