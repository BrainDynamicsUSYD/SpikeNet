function [ OutData ] = ReadH5( files )
% files are cell of array of char array (strings) defining the path of input files
% If given no argument, it searches for matches under CURRENT directory
disp('ReadH5...');
tic;


% Prepare filenames
if nargin == 0
    % If given no argument, search for matches under CURRENT directory
    dir_strut = dir('*out.h5');
    num_files = length(dir_strut);
    files = cell(1,num_files);
    for id_out = 1:num_files
        files{id_out} = dir_strut(id_out).name;
    end
    
else
    % all the corresponding input files
    files_dir = fileparts(files{1}); % assuming all the input files are under the same directory
    if ~isempty(files_dir) % if not under current directory
        in_all_dir = dir(strcat(files_dir,'/*out.h5'));
        in_all = cell(1,1);
        for id_in = 1:length(in_all_dir)
            in_all{id_in} = strcat(files_dir,'/',in_all_dir(id_in).name);
        end
    else % if under current directory
        in_all_dir = dir('/*out.h5');
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
        
                

        % prepare containers
        OutData{id_out}.step_killed = [];
        OutData{id_out}.Num_pop = [];
        OutData{id_out}.dt = [];
        OutData{id_out}.step_tot = [];
        OutData{id_out}.N = [];
        OutData{id_out}.spike_hist = cell(0,0);
        OutData{id_out}.num_spikes = cell(0,0);
        OutData{id_out}.num_ref = cell(0,0);
        OutData{id_out}.ExplVar = [];
        OutData{id_out}.PopPara = cell(0,0);
        OutData{id_out}.SynPara = cell(0,0);
        
        OutData{id_out}.pop_stats.V_mean = cell(0,0);
        OutData{id_out}.pop_stats.V_std = cell(0,0);
        OutData{id_out}.pop_stats.I_input_mean = cell(0,0);
        OutData{id_out}.pop_stats.I_input_std = cell(0,0);
        OutData{id_out}.neuron_stats.IE_ratio = cell(0,0);
        OutData{id_out}.LFP.LFP_neurons = cell(0,0);
        OutData{id_out}.syn_stats = cell(0,0);
        
        
        config_filename = h5read(files{id_out}, '/config_filename/config_filename');
        config_filename = config_filename{1};
        
        
        ev_str = try_h5read(config_filename,'/config/explanatory_variables');
        if ~isempty(ev_str)
            scan_temp = textscan(ev_str{1}, '%s', 'Delimiter', ',');
            scan_temp =  scan_temp{1};
            for s = 1:(length(scan_temp)-1)/2
                str_tmp = strcat('OutData{id_out}.ExplVar.', scan_temp{s*2-1}, '=', scan_temp{s*2},';');
                try
                    eval(str_tmp);
                catch
                end
            end
        end
        
        OutData{id_out}.N = try_h5read(config_filename, '/config/Net/INIT001/N');
        OutData{id_out}.dt = try_h5read(config_filename, '/config/Net/INIT002/dt');
        OutData{id_out}.step_tot = try_h5read(config_filename, '/config/Net/INIT002/step_tot');
        OutData{id_out}.Num_pop = length(OutData{id_out}.N);
        
        OutData{id_out}.step_killed = double(try_h5read(files{id_out}, '/run_away_killed/step')) + 1;
        
        % need debugging!
        for pop_ind = 1:OutData{id_out}.Num_pop
            OutData{id_out}.neuron_sample.neuron_ind{pop_ind,1} = transpose(try_h5read(config_filename, ['/config/pops/pop',num2str(pop_ind-1),'/SAMP001/neurons']));
            OutData{id_out}.neuron_sample.t_ind{pop_ind,1} = transpose(try_h5read(config_filename, ['/config/pops/pop',num2str(pop_ind-1),'/SAMP001/time_points']));
        end

        
        % population results
        for pop_ind = 1:OutData{id_out}.Num_pop
            %
            OutData{id_out}.spike_hist{pop_ind, 1} = transpose(double(try_h5read(files{id_out}, ['/pop_result_' ,num2str(pop_ind-1), '/spike_hist_tot']))) + 1; % be careful here
            OutData{id_out}.num_spikes{pop_ind, 1} = transpose(double(try_h5read(files{id_out}, ['/pop_result_' ,num2str(pop_ind-1), '/num_spikes_pop'])));
            OutData{id_out}.num_ref{pop_ind, 1} = transpose(double(try_h5read(files{id_out}, ['/pop_result_' ,num2str(pop_ind-1), '/num_ref_pop'])));
            %
            OutData{id_out}.PopPara{pop_ind,1} = try_h5read(files{id_out}, ['/pop_result_' ,num2str(pop_ind-1), '/pop_para']);
            %
            OutData{id_out}.pop_stats.V_mean{pop_ind, 1} = transpose(try_h5read(files{id_out}, ['/pop_result_' ,num2str(pop_ind-1), '/stats_V_mean']));
            OutData{id_out}.pop_stats.V_std{pop_ind, 1} = transpose(try_h5read(files{id_out}, ['/pop_result_' ,num2str(pop_ind-1), '/stats_V_std']));
            OutData{id_out}.pop_stats.I_input_mean{pop_ind, 1} = transpose(try_h5read(files{id_out}, ['/pop_result_' ,num2str(pop_ind-1), '/stats_I_input_mean']));
            OutData{id_out}.pop_stats.I_input_std{pop_ind, 1} = transpose(try_h5read(files{id_out}, ['/pop_result_' ,num2str(pop_ind-1), '/stats_I_input_std']));
            OutData{id_out}.pop_stats.I_AMPA_time_avg{pop_ind, 1} = transpose(try_h5read(files{id_out}, ['/pop_result_' ,num2str(pop_ind-1), '/stats_I_AMPA_time_avg']));
            OutData{id_out}.pop_stats.I_NMDA_time_avg{pop_ind, 1} = transpose(try_h5read(files{id_out}, ['/pop_result_' ,num2str(pop_ind-1), '/stats_I_NMDA_time_avg']));
            OutData{id_out}.pop_stats.I_GABA_time_avg{pop_ind, 1} = transpose(try_h5read(files{id_out}, ['/pop_result_' ,num2str(pop_ind-1), '/stats_I_GABA_time_avg']));
            OutData{id_out}.pop_stats.I_ext_time_avg{pop_ind, 1} = transpose(try_h5read(files{id_out}, ['/pop_result_' ,num2str(pop_ind-1), '/stats_I_ext_time_avg']));
            OutData{id_out}.neuron_stats.I_tot_time_mean{pop_ind, 1} = transpose(try_h5read(files{id_out}, ['/pop_result_' ,num2str(pop_ind-1), '/stats_I_tot_time_mean']));
            OutData{id_out}.neuron_stats.I_tot_time_var{pop_ind, 1} = transpose(try_h5read(files{id_out}, ['/pop_result_' ,num2str(pop_ind-1), '/stats_I_tot_time_var']));
            OutData{id_out}.neuron_stats.V_time_mean{pop_ind, 1} = transpose(try_h5read(files{id_out}, ['/pop_result_' ,num2str(pop_ind-1), '/stats_V_time_mean']));
            OutData{id_out}.neuron_stats.V_time_cov{pop_ind, 1} = transpose(try_h5read(files{id_out}, ['/pop_result_' ,num2str(pop_ind-1), '/stats_V_time_cov']));
            OutData{id_out}.neuron_stats.V_time_var{pop_ind, 1} = transpose(try_h5read(files{id_out}, ['/pop_result_' ,num2str(pop_ind-1), '/stats_V_time_var']));
            %
            OutData{id_out}.neuron_stats.IE_ratio{pop_ind, 1} = transpose(try_h5read(files{id_out}, ['/pop_result_' ,num2str(pop_ind-1), '/stats_IE_ratio']));
            %
            OutData{id_out}.LFP.LFP{pop_ind,1} = transpose(try_h5read(files{id_out}, ['/pop_result_' ,num2str(pop_ind-1), '/LFP_data']));
            OutData{id_out}.LFP.LFP_neurons{pop_ind, 1} = transpose(try_h5read(config_filename, ['/config/pops/pop',num2str(pop_ind-1), '/SAMP005/LFP_neurons']));
            %l_tmp = length(OutData{id_out}.LFP.LFP_neurons{pop_ind, 1});
            %N_tmp = OutData{id_out}.N(pop_ind);
            % OutData{id_out}.LFP.LFP_neurons{pop_ind, 1} =  reshape(OutData{id_out}.LFP.LFP_neurons{pop_ind, 1} , [N_tmp l_tmp/N_tmp])';
            
            
        end
        
        
        n_syns = try_h5read(config_filename, '/config/syns/n_syns');
        
        
        if ~isempty(n_syns)
            for syn_ind = 1:n_syns
                %
                OutData{id_out}.syn_stats{syn_ind, 1}.I_mean =  transpose(try_h5read(files{id_out}, ['/syn_result_' ,num2str(syn_ind-1), '/stats_I_mean']));
                OutData{id_out}.syn_stats{syn_ind, 1}.I_std =  transpose(try_h5read(files{id_out}, ['/syn_result_' ,num2str(syn_ind-1), '/stats_std']));
                OutData{id_out}.syn_stats{syn_ind, 1}.s_mean =  transpose(try_h5read(files{id_out}, ['/syn_result_' ,num2str(syn_ind-1), '/stats_s_time_mean']));
                OutData{id_out}.syn_stats{syn_ind, 1}.s_cov =  transpose(try_h5read(files{id_out}, ['/syn_result_' ,num2str(syn_ind-1), '/stats_s_time_cov']));
                OutData{id_out}.syn_stats{syn_ind, 1}.I_time_mean =  transpose(try_h5read(files{id_out}, ['/syn_result_' ,num2str(syn_ind-1), '/stats_I_time_mean']));
                OutData{id_out}.syn_stats{syn_ind, 1}.I_time_var =  transpose(try_h5read(files{id_out}, ['/syn_result_' ,num2str(syn_ind-1), '/stats_I_time_var']));
                
                %
                OutData{id_out}.SynPara{syn_ind, 1} = try_h5read(files{id_out}, ['/syn_result_' ,num2str(syn_ind-1), '/syn_para']);
            end
        end
        
        % neuron sample file
        stamps = file_name(1:end-3);
        for pop_ind = 1:OutData{id_out}.Num_pop
            samp_file = [stamps num2str(pop_ind-1) '_neurosamp.h5'];
            samp_file_mat = [samp_file(1:end-2) 'mat'];
            if exist(samp_file,'file') == 2 && exist(samp_file_mat,'file') ~= 2 % 2 for .mat file
                fprintf('   Generating %s...', samp_file_mat);
                I_AMPA = try_h5read( samp_file,  '/I_AMPA' );
                I_GABA = try_h5read( samp_file,  '/I_GABA' );
                I_K = try_h5read( samp_file,  '/I_K' );
                I_ext = try_h5read( samp_file,  '/I_ext' );
                V = try_h5read( samp_file,  '/V' );
                I_leak = try_h5read( samp_file,  '/I_leak' );
                save(samp_file_mat, 'I_AMPA', 'I_GABA','I_K','I_ext' ,'V', 'I_leak','-v7.3');
                clear 'I_AMPA'  'I_GABA' 'I_K' 'I_ext' 'V''I_leak';
                fprintf('done\n');
            end
        end
        
        
    end
    %                 elseif strfind(tline,'POPD004')
    %                     tline = fgetl(FID);
    %                     scan_temp = textscan(tline,'%f %f','Delimiter',',');
    %                     pop_ind = scan_temp{1}+1; % be careful here!
    %                     sample_size = scan_temp{2};
    %                     tline = fgetl(FID);
    %                     scan_temp = textscan(tline,'%s','Delimiter',',');
    %                     data_name = scan_temp{1};
    %                     for n = 1:length(data_name)
    %                         OutData{id_out}.neuron_sample.(data_name{n}){pop_ind,1} = cell(sample_size,1);
    %                         for sample_ind = 1:sample_size
    %                             tline = fgetl(FID); % read next line
    %                             scan_temp = textscan(tline, '%f', 'Delimiter', ',');
    %                             OutData{id_out}.neuron_sample.(data_name{n}){pop_ind, 1}{sample_ind}= transpose(scan_temp{1});
    %                         end
    %                     end
    %                     for n = 1:length(data_name)
    %                         OutData{id_out}.neuron_sample.(data_name{n}){pop_ind, 1} = cell2mat(OutData{id_out}.neuron_sample.(data_name{n}){pop_ind, 1});
    %                     end
    %                 elseif strfind(tline,'POPD006')
    %                     tline = fgetl(FID);
    %                     scan_temp = textscan(tline,'%f %f %f','Delimiter',',');
    %                     pop_ind = scan_temp{1}+1; % be careful here!
    %                     n_neuron = scan_temp{2};
    %                     n_steps = scan_temp{3};
    %                     tline = fgetl(FID);
    %                     scan_temp = textscan(tline,'%s','Delimiter',',');
    %                     data_name = scan_temp{1};
    %                     data_tmp = zeros(n_neuron*n_steps, length(data_name));
    %                     for i_line = 1:n_neuron*n_steps
    %                         tline = fgetl(FID); % read next line
    %                         scan_temp = textscan(tline, '%f', 'Delimiter', ',');
    %                         data_tmp(i_line,:) = transpose(scan_temp{1});
    %                     end
    %                     for n = 1:length(data_name)
    %                         OutData{id_out}.neuron_sample.(data_name{n}){pop_ind, 1} = transpose(vec2mat(data_tmp(:,n), n_neuron));
    %                     end
    %
    %

    %
    
    %
    %
    %
    %                 elseif strfind(tline,'SYND002')
    %                     if ~isfield(OutData{id_out}, 'syn_sample')
    %                         OutData{id_out}.syn_sample = cell(0,0);
    %                     end
    %                     tline = fgetl(FID);
    %                     scan_temp = textscan(tline, '%d', 'Delimiter', ',');
    %                     OutData{id_out}.syn_sample{end+1,1}.pop_ind_pre = scan_temp{1}(1)+1; % Be careful here! C/C++ index convection!
    %                     OutData{id_out}.syn_sample{end,1}.pop_ind_post = scan_temp{1}(2)+1; % Be careful here! C/C++ index convection!
    %                     OutData{id_out}.syn_sample{end,1}.syn_type = scan_temp{1}(3)+1; % Be careful here! C/C++ index convection!
    %                     sample_size = scan_temp{1}(4);
    %                     for sample_ind = 1:sample_size
    %                         tline = fgetl(FID); % read next line
    %                         scan_temp = textscan(tline, '%f', 'Delimiter', ',');
    %                         OutData{id_out}.syn_sample{end,1}.I(sample_ind,:) = transpose(scan_temp{1});
    %                     end
    %                 elseif strfind(tline,'SAMP002')
    %                     tline = fgetl(FID);
    %                     scan_temp = textscan(tline,'%d','Delimiter',',');
    %                     pop_ind_pre = scan_temp{1}(1)+1; % Be careful here! C/C++ index convection!
    %                     pop_ind_post = scan_temp{1}(2)+1; % Be careful here! C/C++ index convection!
    %                     syn_type = scan_temp{1}(3)+1; % Be careful here! C/C++ index convection!
    %                     for syn_ind = 1:length(OutData{id_out}.syn_sample)
    %                         if OutData{id_out}.syn_sample{syn_ind, 1}.pop_ind_pre == pop_ind_pre && ...
    %                                 OutData{id_out}.syn_sample{syn_ind, 1}.pop_ind_post == pop_ind_post && ...
    %                                 OutData{id_out}.syn_sample{syn_ind, 1}.syn_type == syn_type
    %                             tline = fgetl(FID);
    %                             scan_temp = textscan(tline,'%f','Delimiter',',');
    %                             OutData{id_out}.syn_sample{syn_ind, 1}.neuron_ind = transpose(scan_temp{1} + 1); % Be careful here! C/C++ index convection!
    %                             tline = fgetl(FID);
    %                             scan_temp = textscan(tline,'%f','Delimiter',',');
    %                             OutData{id_out}.syn_sample{syn_ind, 1}.t_ind = find(scan_temp{1}); % extract t_ind for pop_sample
    %                         end
    %                     end
    
    %
    
    %
    %                 elseif strfind(tline,'SYND004')
    %                     if ~isfield(OutData{id_out}, 'syn_tmp_data')
    %                         OutData{id_out}.syn_tmp_data = cell(0,0);
    %                     end
    %                     tline = fgetl(FID);
    %                     scan_temp = textscan(tline, '%d', 'Delimiter', ',');
    %                     OutData{id_out}.syn_tmp_data{end+1,1}.pop_ind_pre = scan_temp{1}(1)+1; % Be careful here! C/C++ index convection!
    %                     OutData{id_out}.syn_tmp_data{end,1}.pop_ind_post = scan_temp{1}(2)+1; % Be careful here! C/C++ index convection!
    %                     OutData{id_out}.syn_tmp_data{end,1}.syn_type = scan_temp{1}(3)+1; % Be careful here! C/C++ index convection!
    %                     data_size = scan_temp{1}(4);
    %                     for sample_ind = 1:data_size
    %                         tline = fgetl(FID); % read next line
    %                         scan_temp = textscan(tline, '%f', 'Delimiter', ',');
    %                         OutData{id_out}.syn_tmp_data{end}.tmp_data(sample_ind,:) = transpose(scan_temp{1});
    %                     end
    %
    
    %                 elseif strfind(tline,'SAMF001')
    %                     if  ~isfield(OutData{id_out},'samp_file')
    %                         OutData{id_out}.samp_file = cell(0,0);
    %                     end
    %                     tline = fgetl(FID);
    %                     scan_temp = textscan(tline, '%s', 'Delimiter', ',');
    %                     OutData{id_out}.samp_file{length(OutData{id_out}.samp_file)+1} = scan_temp{1}{1};
    %
    %                 elseif strfind(tline, 'explanatory variable')
    %                     tline = fgetl(FID);
    %                     scan_temp1 = textscan(tline, '%s', 'Delimiter', ',');
    %                     if strfind(tline,'comment') % read as string
    %                         tline = fgetl(FID);
    %                         scan_temp2 = textscan(tline, '%s', 'Delimiter', ',');
    %                         %eval::OutData{id_out}.ExplVar(end+1,1) =  scan_temp1{1};
    %                         eval(cell2mat(strcat('OutData{id_out}.ExplVar.', scan_temp1{1}, '= cell2mat(scan_temp2{1});')));
    %
    %                     else % read as numeric
    %                         tline = fgetl(FID);
    %                         scan_temp2 = textscan(tline, '%f', 'Delimiter', ',');
    %                         %eval::OutData{id_out}.ExplVar(end+1,1) =  scan_temp1{1};
    %                         eval(cell2mat(strcat('OutData{id_out}.ExplVar.', scan_temp1{1}, '= scan_temp2{1};')));
    %
    %                     end
    %
    %                 elseif strfind(tline, 'INIT004')
    %                 elseif strfind(tline, 'INIT008')
    %                 elseif strfind(tline, 'SAMP005')
    %                 elseif strfind(tline, 'KILL001')
    %                 elseif strfind(tline, 'PARA001')
    %                 elseif strfind(tline, 'INIT005')
    %                 elseif strfind(tline, 'PARA002')
    %                 elseif strfind(tline, 'INIT003')
    %                 elseif strfind(tline, 'SAMP003')
    %                 elseif strfind(tline, 'SAMP004')
    %                 elseif strfind(tline, 'INIT009')
    %                 elseif strfind(tline, 'INIT010')
    %                 elseif strfind(tline, 'INIT011')
    %                 elseif strfind(tline, 'INIT012')
    %                 elseif strfind(tline, 'INIT013')
    %                 elseif strfind(tline, '############')
    %                 elseif strfind(tline, 'MATLAB script')
    %
    %                 else
    %                     warning('unrecognized data type: %s\n', tline);
    %                 end
    %
    
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








