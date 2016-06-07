function ReadYGSamp( files )


% Prepare filenames
if nargin == 0
    % If given no argument, search for matches under CURRENT directory
    dir_strut = dir('*.ygout_samp');
    num_files = length(dir_strut);
    files = cell(1,num_files);
    for id_out = 1:num_files
        files{id_out} = dir_strut(id_out).name;
    end
end

if ~isempty(files)
    % Read ygout file(s) specified by "name" and write data to "OutData"
    for id_out = 1:length(files)
        % OutData{id_out}.file = files{id_out};
        fprintf('Current ReadYGSamp file is: %s\n', files{id_out});
        FID = fopen(files{id_out},'r');
        % prepare containers
        neuron_sample = [];

        while ~feof(FID)
            tline = fgetl(FID);
            % search for data-info line
            if isempty(tline)
                continue;
            elseif strcmp(tline(1), '>')
                
                if strfind(tline,'POPD006')
                    tline = fgetl(FID);
                    scan_temp = textscan(tline,'%f %f %f','Delimiter',',');
                    pop_ind = scan_temp{1}+1; % be careful here!
                    n_neuron = scan_temp{2};
                    n_steps = scan_temp{3};
                    tline = fgetl(FID);
                    scan_temp = textscan(tline,'%s','Delimiter',',');
                    data_name = scan_temp{1};
                    data_tmp = zeros(n_neuron*n_steps, length(data_name));
                    for i_line = 1:n_neuron*n_steps
                        tline = fgetl(FID); % read next line
                        scan_temp = textscan(tline, '%f', 'Delimiter', ',');
                        data_tmp(i_line,:) = transpose(scan_temp{1});
                    end
                    for n = 1:length(data_name)
                        neuron_sample.(data_name{n}){pop_ind} = transpose(vec2mat(data_tmp(:,n), n_neuron));
                    end
                    
                else
                    warning('unrecognized data type in samp file: %s\n', tline);
                end

            end
            
            [~, file_name, ~] = fileparts(files{id_out});
            save([file_name,'_samp.mat'],'neuron_sample')
        end
        fclose(FID);
        
        
    end
    
end % if ~isempty(name)












