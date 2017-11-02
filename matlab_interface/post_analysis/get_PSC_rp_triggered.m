function get_PSC_rp_triggered(stdin)
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

sw_be = 500;
sw_af = 500;

V_rev_I = -80;
V_rev_E = 0;
for f = 1:num_files
    
    samp_file = [files{f}(1:end-11) '0_neurosamp.mat'];
    
    du_af = 1000; % steps
    du_be = 1000;
    
    
    fprintf('\t Loading data from file %s...\n', files{f});
    load(samp_file, 'I_AMPA','I_GABA', 'I_ext', 'V');
    R = load(files{f}, 'LFP', 'step_tot');
    n_neuron = length(R.LFP.wavelet.peak.rp_raw_amp_step );
    EPSC_triggered = cell(n_neuron,1);
    IPSC_triggered = cell(n_neuron,1);
    Vm_triggered = cell(n_neuron,1);
    PSC_rp_amp = cell(n_neuron,1);
    for i = 1:n_neuron
        n_events = length(R.LFP.wavelet.peak.rp_raw_amp_step{i});
        EPSC_triggered{i} = [];
        IPSC_triggered{i} = [];
        Vm_triggered{i} = [];
        PSC_rp_amp{i} = [];
        for j = 1:n_events
            ia = R.LFP.wavelet.peak.rp_raw_amp_step{i}(j); % ripple peak step
            
            % find sharp wave peak step
            if   ia+sw_af-1 <= R.step_tot && ia-sw_be + 1 >=1
                [~,i_sw] = max(R.LFP.LFP_sharpwave(i,ia-sw_be + 1:ia+sw_af-1));
                ia = i_sw + ia-sw_be + 1 - 1;
                
                if   ia+du_af-1 <= R.step_tot && ia-du_be + 1 >=1
                    % patch clamp model
                    EPSC_tmp = (I_AMPA(i,ia-du_be + 1:ia+du_af-1)+I_ext(i,ia-du_be + 1:ia+du_af-1))./(V(i,ia-du_be + 1:ia+du_af-1)-V_rev_E)*(V_rev_I - V_rev_E);
                    EPSC_triggered{i} = [EPSC_triggered{i}; EPSC_tmp]; 
                    IPSC_tmp = I_GABA(i,ia-du_be + 1:ia+du_af-1)./(V(i,ia-du_be + 1:ia+du_af-1)-V_rev_I)*(V_rev_E - V_rev_I);
                    IPSC_triggered{i} = [IPSC_triggered{i}; IPSC_tmp ]; 
                    PSC_rp_amp{i} = [PSC_rp_amp{i}; R.LFP.wavelet.peak.rp_raw_amp{i}(j) ];
                    Vm_triggered{i} = [Vm_triggered{i}; V(i,ia-du_be + 1:ia+du_af-1) ];
                end
            end
        end
    end
    
    save(samp_file,'PSC_rp_amp','-append');
    save(samp_file,'Vm_triggered','-append');
    save(samp_file,'EPSC_triggered', 'IPSC_triggered','-append');
    
end

