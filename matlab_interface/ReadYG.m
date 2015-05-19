function [ OutData ] = ReadYG( files )
% files are cell of array of char array (strings) defining the path of input files
% If given no argument, it searches for matches under CURRENT directory
disp('ReadYG...');
tic;


% Prepare filenames
if nargin == 0
    % If given no argument, search for matches under CURRENT directory
    dir_strut = dir('*.ygout');
    num_files = length(dir_strut);
    files = cell(1,num_files);
    for id_out = 1:num_files
        files{id_out} = dir_strut(id_out).name;
    end
    % all the corresponding input files under CURRENT directory
    in_all_dir = dir('*.ygin');
    in_all = cell(1,1);
    for id_in = 1:length(in_all_dir)
        in_all{id_in} = in_all_dir(id_in).name;
    end
    
else
    % all the corresponding input files
    files_dir = fileparts(files{1}); % assuming all the input files are under the same directory
    if ~isempty(files_dir) % if not under current directory
        in_all_dir = dir(strcat(files_dir,'/*.ygin'));
        in_all = cell(1,1);
        for id_in = 1:length(in_all_dir)
            in_all{id_in} = strcat(files_dir,'/',in_all_dir(id_in).name);
        end
    else % if under current directory
        in_all_dir = dir('*.ygin');
        in_all = cell(1,1);
        for id_in = 1:length(in_all_dir)
            in_all{id_in} = in_all_dir(id_in).name;
        end
    end
end

OutData = cell(1,length(files));

if ~isempty(files)
% Read ygout file(s) specified by "name" and write data to "OutData"
for id_out = 1:length(files)
    % OutData{id_out}.file = files{id_out};
    [file_dir, file_name, ~] = fileparts(files{id_out});
    if ~isempty(file_dir)
        OutData{id_out}.stamp = strcat(file_dir,'/', file_name);
    else
        OutData{id_out}.stamp = file_name;
    end
    fprintf('Current ReadYG file is: %s\n', files{id_out});
    FID = fopen(files{id_out},'r');
    % prepare containers
    OutData{id_out}.step_killed = [];
    OutData{id_out}.Num_pop = [];
    OutData{id_out}.dt = [];
    OutData{id_out}.step_tot = [];
    OutData{id_out}.N = [];
    OutData{id_out}.spike_hist = cell(0,0);
    OutData{id_out}.num_spikes = cell(0,0);
    OutData{id_out}.num_ref = cell(0,0);
    OutData{id_out}.neuron_sample = [];
    OutData{id_out}.syn_sample = cell(0,0);
    OutData{id_out}.syn_stats = cell(0,0);
    OutData{id_out}.ExplVar = [];
    OutData{id_out}.PopPara = cell(0,0);
    OutData{id_out}.SynPara = cell(0,0);
    
    while ~feof(FID)
        tline = fgetl(FID);
        % search for data-info line
        if isempty(tline)
            continue;
        elseif strcmp(tline(1), '>')
            
            
            if strfind(tline,'KILL002');
                tline = fgetl(FID);
                scan_temp = textscan(tline,'%f','Delimiter',','); % using %d will give step_killed a type int32, which will bring trouble!
                OutData{id_out}.step_killed = scan_temp{1} + 1; % Be careful here! C/C++ index convection!
                
                
                
                
                
            elseif strfind(tline,'POPD001')
                tline = fgetl(FID);
                scan_temp = textscan(tline,'%f','Delimiter',',');
                pop_ind = scan_temp{1}+1; % be careful here!
                tline = fgetl(FID);
                scan_temp = textscan(tline, '%f', 'Delimiter', ',');
                OutData{id_out}.spike_hist{pop_ind,1} = transpose(scan_temp{1} + 1); % Be careful here! C/C++ index convection!
                tline = fgetl(FID);
                scan_temp = textscan(tline, '%f', 'Delimiter', ',');
                OutData{id_out}.num_spikes{pop_ind,1} = transpose(scan_temp{1});
                tline = fgetl(FID);
                scan_temp = textscan(tline, '%f', 'Delimiter', ',');
                OutData{id_out}.num_ref{pop_ind,1} = transpose(scan_temp{1});
                
                
                
             elseif strfind(tline,'POPD003')
                tline = fgetl(FID);
                scan_temp = textscan(tline,'%d','Delimiter',',');
                pop_ind = scan_temp{1}+1; % be careful here!
                tline = fgetl(FID);
                scan_temp = textscan(tline,'%f','Delimiter',',');
                V_mean = transpose(scan_temp{1});
                tline = fgetl(FID);
                scan_temp = textscan(tline,'%f','Delimiter',',');
                V_std = transpose(scan_temp{1});
                tline = fgetl(FID);
                scan_temp = textscan(tline,'%f','Delimiter',',');
                I_input_mean = transpose(scan_temp{1});
                tline = fgetl(FID);
                scan_temp = textscan(tline,'%f','Delimiter',',');
                I_input_std = transpose(scan_temp{1});
                OutData{id_out}.pop_stats.V_mean{pop_ind} = V_mean;
                OutData{id_out}.pop_stats.V_std{pop_ind} = V_std;
                OutData{id_out}.pop_stats.I_input_mean{pop_ind} = I_input_mean;
                OutData{id_out}.pop_stats.I_input_std{pop_ind} = I_input_std;
                
                
                       
            elseif strfind(tline,'POPD004')
                tline = fgetl(FID);
                scan_temp = textscan(tline,'%f %f','Delimiter',',');
                pop_ind = scan_temp{1}+1; % be careful here!
                sample_size = scan_temp{2};
                tline = fgetl(FID);
                scan_temp = textscan(tline,'%s','Delimiter',',');
                data_name = scan_temp{1};
                for n = 1:length(data_name)
                    OutData{id_out}.neuron_sample.(data_name{n}){pop_ind,1} = cell(sample_size,1);
                    for sample_ind = 1:sample_size
                        tline = fgetl(FID); % read next line
                        scan_temp = textscan(tline, '%f', 'Delimiter', ',');
                        OutData{id_out}.neuron_sample.(data_name{n}){pop_ind}{sample_ind}= transpose(scan_temp{1});
                    end
                end
                for n = 1:length(data_name)
                    OutData{id_out}.neuron_sample.(data_name{n}){pop_ind} = cell2mat(OutData{id_out}.neuron_sample.(data_name{n}){pop_ind});
                end
            elseif strfind(tline,'SAMP001')
                tline = fgetl(FID);
                scan_temp = textscan(tline,'%f','Delimiter',',');
                pop_ind = scan_temp{1}(1)+1; % be careful here!
                fgetl(FID); % skip line: data type
                tline = fgetl(FID);
                scan_temp = textscan(tline,'%f','Delimiter',',');
                OutData{id_out}.neuron_sample.neuron_ind{pop_ind,1} = transpose(scan_temp{1} + 1); % Be careful here! C/C++ index convection!
                tline = fgetl(FID);
                scan_temp = textscan(tline,'%f','Delimiter',',');
                OutData{id_out}.neuron_sample.t_ind{pop_ind,1} = find(scan_temp{1}); % extract t_ind for pop_sample
                
                
                
               
            elseif strfind(tline,'SYND002')
                tline = fgetl(FID);
                scan_temp = textscan(tline, '%d', 'Delimiter', ',');
                OutData{id_out}.syn_sample{end+1,1}.pop_ind_pre = scan_temp{1}(1)+1; % Be careful here! C/C++ index convection!
                OutData{id_out}.syn_sample{end,1}.pop_ind_post = scan_temp{1}(2)+1; % Be careful here! C/C++ index convection!
                OutData{id_out}.syn_sample{end,1}.syn_type = scan_temp{1}(3)+1; % Be careful here! C/C++ index convection!
                sample_size = scan_temp{1}(4);
                for sample_ind = 1:sample_size
                        tline = fgetl(FID); % read next line
                        scan_temp = textscan(tline, '%f', 'Delimiter', ',');
                        OutData{id_out}.syn_sample{end,1}.I(sample_ind,:) = transpose(scan_temp{1});
                end
            elseif strfind(tline,'SAMP002')
                tline = fgetl(FID);
                scan_temp = textscan(tline,'%d','Delimiter',',');
                pop_ind_pre = scan_temp{1}(1)+1; % Be careful here! C/C++ index convection!
                pop_ind_post = scan_temp{1}(2)+1; % Be careful here! C/C++ index convection!
                syn_type = scan_temp{1}(3)+1; % Be careful here! C/C++ index convection!
                for syn_ind = 1:length(OutData{id_out}.syn_sample)
                    if OutData{id_out}.syn_sample{syn_ind}.pop_ind_pre == pop_ind_pre && ...
                            OutData{id_out}.syn_sample{syn_ind}.pop_ind_post == pop_ind_post && ...
                            OutData{id_out}.syn_sample{syn_ind}.syn_type == syn_type
                        tline = fgetl(FID);
                        scan_temp = textscan(tline,'%f','Delimiter',',');
                        OutData{id_out}.syn_sample{syn_ind}.neuron_ind = transpose(scan_temp{1} + 1); % Be careful here! C/C++ index convection!
                        tline = fgetl(FID);
                        scan_temp = textscan(tline,'%f','Delimiter',',');
                        OutData{id_out}.syn_sample{syn_ind}.t_ind = find(scan_temp{1}); % extract t_ind for pop_sample
                    end
                end

                
            elseif strfind(tline,'SYND003')
                tline = fgetl(FID);
                scan_temp = textscan(tline, '%d', 'Delimiter', ',');
                OutData{id_out}.syn_stats{end+1,1}.pop_ind_pre = scan_temp{1}(1)+1; % Be careful here! C/C++ index convection!
                OutData{id_out}.syn_stats{end,1}.pop_ind_post = scan_temp{1}(2)+1; % Be careful here! C/C++ index convection!
                OutData{id_out}.syn_stats{end,1}.syn_type = scan_temp{1}(3)+1; % Be careful here! C/C++ index convection!
                tline = fgetl(FID);
                scan_temp = textscan(tline,'%f','Delimiter',',');
                I_mean = transpose(scan_temp{1});
                tline = fgetl(FID);
                scan_temp = textscan(tline,'%f','Delimiter',',');
                I_std = transpose(scan_temp{1});
                OutData{id_out}.syn_stats{end}.I_mean = I_mean;
                OutData{id_out}.syn_stats{end}.I_std = I_std;
                
                
            elseif strfind(tline,'POPD002')
                tline = fgetl(FID);
                scan_temp = textscan(tline,'%f %f','Delimiter',',');
                pop_ind = scan_temp{1}+1; % be careful here!
                num_para = scan_temp{2};
                for p = 1:num_para
                    tline = fgetl(FID); % read next line
                    scan_temp = textscan(tline, '%s %f', 'Delimiter', ',');
                    para_name = scan_temp{1}{1};
                    para_value = scan_temp{2};
                    OutData{id_out}.PopPara{pop_ind,1}.(para_name) = para_value;
                end
            elseif strfind(tline,'SYND001')
                tline = fgetl(FID);
                scan_temp = textscan(tline, '%f %s', 'Delimiter', ',');
                num_para = scan_temp{1}; % number of parameters
                num_syn = length(OutData{id_out}.SynPara)+1;
                for p = 1:num_para
                    tline = fgetl(FID); % read next line
                    scan_temp = textscan(tline, '%s %f', 'Delimiter', ',');
                    para_name = scan_temp{1}{1};
                    para_value = scan_temp{2};
                    OutData{id_out}.SynPara{num_syn,1}.(para_name) = para_value;
                end
                
                
                
                


                
                
            elseif strfind(tline,'INIT002')
                tline = fgetl(FID);
                scan_temp = textscan(tline, '%f %f', 'Delimiter', ',');
                OutData{id_out}.dt = scan_temp{1};
                OutData{id_out}.step_tot = scan_temp{2};
            elseif strfind(tline, 'INIT001')
                tline = fgetl(FID);
                scan_temp = textscan(tline, '%f', 'Delimiter', ',');
                OutData{id_out}.N = scan_temp{1};
                OutData{id_out}.Num_pop = length(OutData{id_out}.N);
                
                
                
                
                
                
                
                
            elseif strfind(tline, 'explanatory variable')
                tline = fgetl(FID);
                scan_temp1 = textscan(tline, '%s', 'Delimiter', ',');
                if strfind(tline,'comment') % read as string
                    tline = fgetl(FID);
                    scan_temp2 = textscan(tline, '%s', 'Delimiter', ',');
                    %eval::OutData{id_out}.ExplVar(end+1,1) =  scan_temp1{1};
                    eval(cell2mat(strcat('OutData{id_out}.ExplVar.', scan_temp1{1}, '= cell2mat(scan_temp2{1});')));
                    
                else % read as numeric
                    tline = fgetl(FID);
                    scan_temp2 = textscan(tline, '%f', 'Delimiter', ',');
                    %eval::OutData{id_out}.ExplVar(end+1,1) =  scan_temp1{1};
                    eval(cell2mat(strcat('OutData{id_out}.ExplVar.', scan_temp1{1}, '= scan_temp2{1};')));
                    
                end
                
                
                
                
                
            else
                warning('unrecognized data type: %s\n', tline);
            end
            
            
        end
    
    end
    fclose(FID);

        
end




% Re-formatting data
fprintf('\t Re-formatting data...\n');
Result_num = length(OutData);
for r_num = 1:Result_num
    
        
    % adaptation of step_killed
    if OutData{r_num}.step_killed > 0
        OutData{r_num}.step_tot = OutData{r_num}.step_killed;
    end
    
%     % Re-format VI_sample into coloum matrix
%     OutData{r_num} = ReformatVI_sample(OutData{r_num});

    % Reformat spike history data
    OutData{r_num} = ReformatSpikeHistory(OutData{r_num});
 
    % Discard transient data
    OutData{r_num} = DiscardTransientData(OutData{r_num});
    
    % Reduce solution
    OutData{r_num} = ReduceSolution(OutData{r_num});
end

end % if ~isempty(name)

toc;

end


% function Data = ReformatVI_sample(Data)
%     % Re-format VI_sample into coloum matrix
%     if any( strcmp(fieldnames(Data), 'VI_sample') ) && ~isempty(Data.VI_sample) 
%         fname_cell = fieldnames(Data.VI_sample);
%         for f = 1:length(fname_cell)
%             fn = fname_cell{f};
%             Data.VI_sample.(fn) = cell2mat(Data.VI_sample.(fn));
%         end
%     end
% end








