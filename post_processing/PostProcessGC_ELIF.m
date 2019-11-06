function PostProcessGC_ELIF(stdin)
% stdin is paths of input files in unix-style, i.e., separated by spaces
% If given no argument, it searches for matches under CURRENT directory


% Prepare files
if nargin == 0
    dir_strut = dir('*out.h5');
    num_files = length(dir_strut);
    files = cell(1,num_files);
    for id_out = 1:num_files
        files{id_out} = dir_strut(id_out).name;
    end
else
    % stdin, i.e., file pathes and names separated by space
    sfiles = textscan(stdin,'%s'); % cell array of file path+names
    num_files = size(sfiles{1},1);
    for i = 1:num_files
        files{i} = cell2mat(sfiles{1}(i,:));
    end
end

% save figures
save_fig = 1; % -1 for no figure, 0 for displaying figure, 1 for saving figure
% Start processing
for id_out = 1:num_files
    % start from .ygout  files
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    R = ReadH5( files(id_out) ); % read .ygout file into matlab data struct
    R = AnalyseGC_ELIF(R); % do some simple analysis
    SaveRYG(R);
    disp('Done');
    %RasterGC_ELIF(R, save_fig); % generate raster plot for spiking history
    %CascadeAnalysis_GC(R)
    %GCplotV_onegrid(R); % record the video of V combing E and I
    %     GCplotV_onegrid_nosquare(R)
    %     GCplotV_onegrid_savetransientfig(R); % save the shot of specific duration.
    %GCtrackPattern(R,1); % track patterns
    % GCplotSpike_onegrid(R); % record the video of spike
    %         GCplotCurr_onegrid(R,'AMPA'); % record the video of excitatory current
    %         GCplotCurr_onegrid(R,'GABA'); % record the video of inhibitory current
    
    get_initialV(R,R{1}.step_tot/50); % get the last V for external setting V
%     try
%         CombineSampleData(R)
%     catch
%     end
    %     GCplot_totalCurr_onegrid(R)
%     get_combinedSpikeHist(R)
end
%     CombineFig(1);% 1: firing rate and raster
end

