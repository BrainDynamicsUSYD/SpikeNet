function get_V_LFP_rp_triggered(stdin)
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


for f = 1:num_files

    samp_file = [files{f}(1:end-11) '0_neurosamp.mat'];
    
    du_af = 4000; % steps
    du_be = 4000;
    V_rp_triggered = [];
    LFP_rp_triggered = [];
    
    fprintf('\t Loading data from file %s...\n', files{f});
    load(samp_file, 'V');
    R = load(files{f}, 'LFP', 'step_tot');
    for i = 1:length(R.LFP.wavelet.peak.rp_raw_amp_step)
        for j = 1:length(R.LFP.wavelet.peak.rp_raw_amp_step{i} )
            ia = R.LFP.wavelet.peak.rp_raw_amp_step{i}(j);
            
            if   ia+du_af-1 <= R.step_tot && ia-du_be + 1 >=1
                V_rp_triggered = [V_rp_triggered; V(i,ia-du_be + 1:ia+du_af-1)]; %#ok<AGROW>
                LFP_rp_triggered = [LFP_rp_triggered; R.LFP.LFP_broad(i,ia-du_be + 1:ia+du_af-1)]; %#ok<AGROW>
            end
        end
    end
    
    
    save(samp_file,'V_rp_triggered', 'LFP_rp_triggered','-append');
end

