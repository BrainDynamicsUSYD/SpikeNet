function PostProcessYG(stdin)
% stdin is paths of input files in unix-style, i.e., separated by spaces
% If given no argument, it searches for matches under CURRENT directory


% Prepare files
if nargin == 0
    dir_strut = dir('*.ygout');
    num_files = length(dir_strut);
    file_type = 1;
    if num_files == 0
        dir_strut = dir('*out.h5');
        num_files = length(dir_strut);
        file_type = 2;
    end
    
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
    [~,~,app] = fileparts(files{1});
    if strcmp(app, '.ygout')
        file_type = 1;
    elseif strcmp(app, '.h5')
        file_type = 2;
    end
end

% save figures
save_fig = 1; % -1 for no figure, 0 for displaying figure, 1 for saving figure
% Start processing
for id_out = 1:num_files
    % start from .ygout  files
    fprintf('Processing output file No.%d out of %d...\n', id_out, num_files);
    fprintf('\t File name: %s\n', files{id_out});
    if  file_type == 1
        R = ReadYG( files(id_out) ); % read .ygout file into matlab data struct
    elseif file_type == 2
        R = ReadH5( files(id_out) ); % read .ygout file into matlab data struct
    end
    R = AnalyseYG(R); % do some simple analysis
    
    try
        R = get_SWR(R{1});R = {R};
    catch
    end
    
    
    SaveRYG(R);
    disp('Done');
    for ind = 1:length(R)
        if isfield(R{ind}, 'samp_file')
            Read_and_save_YGSamp(R{ind}.samp_file, R{ind});
        end
    end
    
    % RasterYG(R, save_fig); % generate raster plot for spiking history

end


end

