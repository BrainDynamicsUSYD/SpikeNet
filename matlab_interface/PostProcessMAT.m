function PostProcessMAT(stdin)
% stdin is paths of input files in unix-style, i.e., separated by spaces
% If given no argument, it searches for matches under CURRENT directory


% Prepare files
if nargin == 0
    dir_strut = dir('*.mat');
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

% save figures
save_fig = 1; % -1 for no figure, 0 for displaying figure, 1 for saving figure
% Start processing
for id_out = 1:num_files
    
    % start form .mat files
    fprintf('Loading RYG.mat file %s...', files{i});
    R_temp = load(files{i});
    disp('done.');
    %%%%%%% do something here
    
%     R_temp = cluster_sorted_rate(R_temp);
%     cluster = R_temp.cluster;
%     save(files{i},'cluster', '-append');
    
%     R_temp = get_CC_pop(R_temp);
%     R_temp = get_EI_current_crosscorr(R_temp);
%     Balance = R_temp.Balance;
%     Analysis = R_temp.Analysis;
%     save(files{i},'Balance', 'Analysis', '-append');

    R_temp = cluster_sorted_rate(R_temp);
    cluster = R_temp.cluster;
    save(files{i},'cluster', '-append');

    
    % save(files{i},'-struct', 'R_temp', '-v7.3'); % -v7.3 for >2GB
   % RasterYG({R_temp}, save_fig);
    
   
   
    %%%%%%% do something above
    disp('Data processed and saved.');
    
end


end
