function Read_and_save_YGSamp( files, R )


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
        FID = fopen(files{id_out},'r');
        % prepare containers
        neuron_sample = [];
        saved = 0;
        if FID == -1
            continue;
        end
        fprintf('Processing ygout_samp file No.%d out of %d...\n', id_out, length(files));
        fprintf('\t Current ReadYGSamp file is: %s\n', files{id_out});
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
                    % if the simulation is killed, correct n_steps
                    if nargin == 2
                        if R.step_killed > 0
                            n_steps = sum(R.neuron_sample.t_ind{pop_ind} <= R.step_killed);
                        end
                    else
                        warning('Need the 2nd input argument to correct n_steps if the simulation is killed.')
                    end
                    
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
            
            [path, file_name, ~] = fileparts(files{id_out});
            save([path, file_name,'_samp.mat'],'neuron_sample')
            saved = 1;
        end
        fclose(FID);
        if saved == 1
            delete(files{id_out});
        end
        %fclose(FID);
        disp('Done');
        
    end
    
end % if ~isempty(name)












