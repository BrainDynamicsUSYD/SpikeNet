function PostProcessMAT(stdin)
% stdin is paths of input files in unix-style, i.e., separated by spaces
% If given no argument, it searches for matches under CURRENT directory


% Prepare files
if nargin == 0
    dir_strut = dir('*.mat');
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
    
    
    
    %     R_temp = rmfield(R_temp,'cluster');
    %     R_temp = cluster_sorted_rate(R_temp);
    %     cluster = R_temp.cluster;
    %     save(files{i},'cluster', '-append');
    
    
    % R_temp = get_CC_network(R_temp);
    % R_temp = get_CC_pop(R_temp, 1);
    % R_temp = get_CV2_ISI(R_temp);
    %     R_temp = get_ISI_low_high(R_temp);
    %     Analysis = R_temp.Analysis;
    %     save(files{i},'Analysis', '-append');
    
    
%     R_temp = avalanche_detect(R_temp);
%     avalanche = R_temp.avalanche;
%     save(files{i},'avalanche', '-append');
    
    
%     R_temp = get_grid_firing_centre(R_temp);
%     grid = R_temp.grid;
%     save(files{i},'grid', '-append');
    
    R_temp = get_stPR(R_temp);
    stPR = R_temp.stPR;
    save(files{i},'stPR', '-append');
    

    % R_temp = rmfield(R_temp,{'C_rate','C_label','up_down_analysis'});
    
    
    %     R_temp = get_CC_pop(R_temp);
    %     R_temp = get_EI_current_crosscorr(R_temp);
    %     Balance = R_temp.Balance;
    %     Analysis = R_temp.Analysis;
    %     save(files{i},'Balance', 'Analysis', '-append');
    
    %     R_temp = cluster_sorted_rate(R_temp);
    %     cluster = R_temp.cluster;
    %     save(files{i},'cluster', '-append');
    %     ClusterRasterYG({R_temp}, save_fig);
    %     RasterYG({R_temp}, save_fig);
    
    %save(files{i},'-struct', 'R_temp', '-v7.3'); % -v7.3 for >2GB
    
    %     R_temp = get_neuron_sample_stats(R_temp);
    %     neuron_sample_stats = R_temp.neuron_sample_stats;
    %     save(files{i},'neuron_sample_stats', '-append');
    
    
    
    %%%%%%% do something above
    disp('Data processed and saved.');
    
end


end
