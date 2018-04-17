function get_g_spike_triggered(stdin)
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

sp_be = 500; % 500 dt = 50 ms
sp_af = 0;

for f = 1:num_files
    
    samp_file = [files{f}(1:end-11) '0_neurosamp.mat'];

    fprintf('\t Loading data from file %s...\n', files{f});
    load(samp_file, 'I_ext', 'I_AMPA','I_GABA', 'V');
    
    gI = -I_GABA./(V+80);
    gE = -(I_AMPA+I_ext)./V;
    
    R = load(files{f}, 'step_tot','neuron_sample','spike_hist');
    ind =  R.neuron_sample.neuron_ind{1};
    n_neuron = length(ind);
    gE_sp_triggered = [];
    gI_sp_triggered = [];
    sh = R.spike_hist{1}(ind,:);
    for i = 1:n_neuron
        n_events = sum(sh(i,:));
        sp_tmp = find(sh(i,:));
        for j = 1:10:n_events
            ia = sp_tmp(j); % ripple peak step
            
            % find sharp wave peak step
            if  ia+sp_af-1 <= R.step_tot && ia-sp_be + 1 >=1
                gE_sp_triggered = [gE_sp_triggered; gE(i,ia-sp_be + 1:10:ia+sp_af-1 )]; %#ok<*AGROW>
                gI_sp_triggered = [gI_sp_triggered; gI(i,ia-sp_be + 1:10:ia+sp_af-1 )];
            end
        end
    end
    
    gE_sp_triggered = mean(gE_sp_triggered);
    gI_sp_triggered = mean(gI_sp_triggered);
    save(samp_file,'gE_sp_triggered', 'gI_sp_triggered','-append');
    
end

