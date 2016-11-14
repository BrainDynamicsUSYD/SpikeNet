function [ R ] = get_SWR( R )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% ref: Hippocampal place-cell sequences depict future paths to remembered goals.

ripple_event.spike_sort_range = 6;
ripple_event.no_std = 1.5; % 2 std above the mean
ripple_event.peak_no_std = 3; % 4 std above the mean
ripple_event.ripple_min_ms = 15;

dt = R.dt;
fs = 1/(dt*1e-3); % sampling frequency (Hz)

try
    LFP = R.LFP.LFP{1};
    [no, steps] = size(LFP);
    % no = 1; % for testing
    
    
    
    
    % Butterworth filter
    order = 4; % 4th order
    lowFreq_br = 1; % broad band (1-1000 Hz)
    hiFreq_br = 1000;
    Wn = [lowFreq_br hiFreq_br]/(fs/2);
    [b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
    
    for i = 1:no
        LFP_broad(i,:) = filter(b,a,LFP(i,:)); %#ok<AGROW>
    end
    
    % Butterworth filter
    order = 4; % 4th order
    lowFreq_sp = 1; % broad band (10?250 Hz)
    hiFreq_sp = 50;
    Wn = [lowFreq_sp hiFreq_sp]/(fs/2);
    [b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
    
    for i = 1:no
        LFP_sharpwave(i,:) = filter(b,a,LFP(i,:)); %#ok<AGROW>
    end
    
    
    % Butterworth filter
    order = 4; % 4th order
    R.LFP.lowFreq = 100; % ripple band (default values for this function are 150-250 Hz)
    R.LFP.hiFreq = 250;
    Wn = [R.LFP.lowFreq R.LFP.hiFreq]/(fs/2);
    [b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
    gaus_width = 5; %ms
    [ Kernel ] = spike_train_kernel_YG( gaus_width, dt, 'gaussian_unit' );
    for i = 1:no
        LFP_ripple(i,:) = filter(b,a,LFP(i,:)); %#ok<AGROW>
        % hilbert transformation & gaussian smoothing
        LFP_ripple_hilbert(i,:) = conv(abs(hilbert(LFP_ripple(i,:))), Kernel,'same'); %#ok<AGROW>
        %     rms_window_ms = 17; %ms
        %     window_steps = round(rms_window_ms/dt);
        %     Kernel_rms = ones(1,window_steps)/window_steps;
        %     LFP_ripple_rms(i,:) = sqrt(conv(LFP_ripple(i,:).^2, Kernel_rms,'same')); %#ok<AGROW>
    end
    
    
    % Output results
    R.LFP.LFP_broad = LFP_broad;
    R.LFP.LFP_ripple = LFP_ripple;
    R.LFP.LFP_sharpwave = LFP_sharpwave;
    R.LFP.LFP_ripple_hilbert = LFP_ripple_hilbert;
    R.LFP.gauss_width = gaus_width;
    clear LFP_broad LFP_ripple  LFP_sharpwave LFP_ripple_hilbert;
catch EM
    EM
end
% % Discard transient data
% transient_ms = 0; %ms;
% R.LFP.transient_steps = transient_steps;
% transient_steps = round(transient_ms/dt);
% if steps <= transient_steps
%     error('Not enough LFP data!')
% end
% LFP_broad = LFP_broad(:,transient_steps+1:end);
% LFP_ripple = LFP_ripple(:,transient_steps+1:end);
% LFP_ripple_hilbert = LFP_ripple_hilbert(:,transient_steps+1:end);


% Ripple event detection
[no, steps] = size(R.LFP.LFP_broad);
spike_hist = R.spike_hist{1};
spike_tot = sum(sum(spike_hist));

ripple_min_steps = round(ripple_event.ripple_min_ms/dt);

% Discard transient data
transient_ms = 200; %ms;
ripple_event.transient_ms = transient_ms;
transient_steps = min(round(transient_ms/dt), steps);
ripple_event.index1 = cell(1,no);
ripple_event.index2 = cell(1,no);
ripple_event.index3 = cell(1,no);
ripple_event.inside_rate = cell(1,no);
ripple_event.outside_rate = cell(1,no);
ripple_event.Hz = zeros(1,no);
ripple_event.is_SWR = zeros(size(R.LFP.LFP_ripple_hilbert));

iter_num = 100;
hil_mean_baseline_hist = zeros(no,iter_num);
hil_std_baseline_hist = zeros(no,iter_num);
for i = 1:no
    LFP_ripple_hilbert_tmp = R.LFP.LFP_ripple_hilbert(i,:);
    is_SWR_tmp = zeros(size(LFP_ripple_hilbert_tmp));
    for r_iter = 1:iter_num
        hil_mean_baseline_hist(i,r_iter) =  mean(LFP_ripple_hilbert_tmp(~is_SWR_tmp));
        hil_std_baseline_hist(i,r_iter)  =  std(LFP_ripple_hilbert_tmp(~is_SWR_tmp));
        
        if r_iter > 1
            mean_diff = abs(hil_mean_baseline_hist(i,r_iter) - hil_mean_baseline_hist(i,r_iter-1))/hil_mean_baseline_hist(i,r_iter-1);
            std_diff = abs(hil_std_baseline_hist(i,r_iter) - hil_std_baseline_hist(i,r_iter-1))/hil_std_baseline_hist(i,r_iter-1);
            if mean_diff < 0.05 && std_diff < 0.05
                ripple_event.is_SWR(i,:) = is_SWR_tmp;
                break
            end
        end
        
        % thresholding
        is_SWR_tmp = LFP_ripple_hilbert_tmp - hil_mean_baseline_hist(i,r_iter) > ripple_event.no_std*hil_std_baseline_hist(i,r_iter);
        
        % Discard transient data
        is_SWR_tmp(1:transient_steps) = 0;
        
        % persistence_requirement
        is_SWR_tmp = persistence_requirement( is_SWR_tmp, ripple_min_steps );
        
        % peak height requirement
        peak_min = hil_mean_baseline_hist(i,r_iter) + ripple_event.peak_no_std*hil_std_baseline_hist(i,r_iter);
        [ is_SWR_tmp ] = peak_height_requirement( is_SWR_tmp, LFP_ripple_hilbert_tmp, peak_min );
        
    end
    
    cutoff = 1;
    [~, ripple_du, flat_du, ripple_start, ~] = seq_postprocess(ripple_event.is_SWR(i,:), 1, cutoff);
    ripple_event.Hz(i) = length(ripple_du)/(R.dt*R.step_tot*10^-3);
    ripple_event.ripple_du_steps{i} = ripple_du;
    ripple_event.flat_du_steps{i} = flat_du;
    ripple_event.ripple_start_steps{i} = ripple_start;
    
    % LFP_neurons = logical(R.LFP.LFP_neurons{1}(i,:));
    s_tmp = R.ExplVar.LFP_range_sigma;
    spike_sort_neurons = R.LFP.LFP_neurons{1}(i,:) >= 1/(s_tmp*sqrt(2*pi))*exp(-0.5*(ripple_event.spike_sort_range/s_tmp)^2);
    
    % firing rate in and out of SWR
    ripple_event.inside_rate{i} = transpose(sum(spike_hist(spike_sort_neurons, is_SWR_tmp), 2)/(R.dt*sum(is_SWR_tmp)*10^-3));
    ripple_event.outside_rate{i} = transpose(sum(spike_hist(spike_sort_neurons, ~is_SWR_tmp), 2)/(R.dt*sum(~is_SWR_tmp)*10^-3));
    
    
    % three indexes for population synchrony
    % ref: Preconfigured, skewed distribution of firing rates in the hippocampus and entorhinal cortex
    N_events = length(ripple_du);
    ripple_event.index1{i} = zeros(1,N_events);
    ripple_event.index2{i} = zeros(1,sum(spike_sort_neurons));
    ripple_event.index3{i} = zeros(1,sum(spike_sort_neurons));
    for r = 1:N_events
        spike_tmp = spike_hist(spike_sort_neurons,ripple_start(r):ripple_start(r)+ripple_du(r)-1);
        ripple_event.index1{i}(r) =  sum(sum(spike_tmp))/spike_tot; % percentage of spikes in each SWR
        participated = (sum(spike_tmp,2))>0;
        ripple_event.index2{i}(participated) = ripple_event.index2{i}(participated) + 1/N_events; % number of neurons participated in each SWR
        ripple_event.index3{i} = ripple_event.index3{i} + (sum(spike_tmp,2)/N_events)'; % mean number of spikes per SWR
    end
end

ripple_event.hil_mean_baseline_hist = hil_mean_baseline_hist;
ripple_event.hil_std_baseline_hist = hil_std_baseline_hist;
ripple_event.hil_mean_baseline = zeros(no,1);
ripple_event.hil_std_baseline = zeros(no,1);
for i = 1:no
    [~,~,sw_amp_mean_base_tmp] = find(ripple_event.hil_mean_baseline_hist(i,:), 1, 'last');
    [~,~,sw_amp_std_base_tmp] = find(ripple_event.hil_std_baseline_hist(i,:), 1, 'last');
    ripple_event.hil_mean_baseline(i) = sw_amp_mean_base_tmp;
    ripple_event.hil_std_baseline(i) = sw_amp_std_base_tmp;
end


%
freqrange = [R.LFP.lowFreq R.LFP.hiFreq];
Fs = 1000/dt;
fc = centfrq('cmor1.5-1');
scalerange = fc./(freqrange*(1/Fs));
scales = scalerange(end):0.5:scalerange(1);
pseudoFreq = scal2frq(scales,'cmor1.5-1',1/Fs); % pseudo-frequencies
% coeffs = cell(1,no);
peak.rp_freq = cell(1,no);
peak.rp_freq_step = cell(1,no);
peak.rp_wl_amp = cell(1,no);
peak.rp_raw_amp = cell(1,no);
peak.rp_raw_amp_step = cell(1,no);
peak.rp_amp_no_std = cell(1,no);
peak.sw_amp = cell(1,no);
peak.prn = cell(1,no);
sw_average = cell(1,no);
rp_average = cell(1,no);
for i = 1:no
    x_tmp = R.LFP.LFP_ripple(i,:);
    coeffs_tmp = abs(cwt(x_tmp,scales,'cmor1.5-1'))';
    % coeffs{i} = coeffs_tmp;
    
    
    peak.sw_amp{i} = [];
    for j =  1:length(ripple_event.ripple_du_steps{i})
        a_tmp = ripple_event.ripple_start_steps{i}(j);
        l_tmp = ripple_event.ripple_du_steps{i}(j);
        
        [max_coffs_tmp, freq_ind_tmp] = max(coeffs_tmp(a_tmp:a_tmp+l_tmp-1, :), [], 2);
        [rp_wl_amp, ind_tmp] = max(max_coffs_tmp);
        peak.rp_freq{i} = [peak.rp_freq{i} pseudoFreq(freq_ind_tmp(ind_tmp))];
        peak.rp_freq_step{i} = [peak.rp_freq_step{i} a_tmp+ind_tmp-1];
        peak.rp_wl_amp{i} = [peak.rp_wl_amp{i} rp_wl_amp];
        
        [rp_raw_amp, rp_raw_amp_step] = max(R.LFP.LFP_ripple(i, a_tmp:a_tmp+l_tmp-1));
        peak.rp_raw_amp{i} = [peak.rp_raw_amp{i} rp_raw_amp];
        peak.rp_raw_amp_step{i} = [peak.rp_raw_amp_step{i} a_tmp+rp_raw_amp_step-1];

    end
    peak.rp_amp_no_std{i} =  (peak.rp_raw_amp{i} - ripple_event.hil_mean_baseline(i)) / ripple_event.hil_std_baseline(i);
      
    % sharp-wave and ripple correlation
    swr_hw = 75; % ms
    swr_hw_steps = round(swr_hw/dt);
    sw_average{i} = zeros(1,swr_hw_steps*2+1);
    rp_average{i} = zeros(1,swr_hw_steps*2+1);
    for j = 1:length(peak.rp_raw_amp_step{i})
        peak_t = peak.rp_raw_amp_step{i}(j);
        if swr_hw_steps + peak_t < R.step_tot && -swr_hw_steps + peak_t > 0
            peak.sw_amp{i} = [peak.sw_amp{i} max(R.LFP.LFP_sharpwave(i, -swr_hw_steps + peak_t: swr_hw_steps + peak_t))];
            peak.prn{i} = [ peak.prn{i} min(R.LFP.LFP_broad(i, peak_t: swr_hw_steps + peak_t))];
            sw_average{i} = sw_average{i} + R.LFP.LFP_sharpwave(i, -swr_hw_steps + peak_t: swr_hw_steps + peak_t);
            rp_average{i} = rp_average{i} + R.LFP.LFP_ripple(i, -swr_hw_steps + peak_t: swr_hw_steps + peak_t);
        else
            peak.sw_amp{i} = [peak.sw_amp{i} NaN];
            peak.prn{i} = [ peak.prn{i} NaN];
        end
    end
    sw_average{i} = sw_average{i}/length(peak.rp_raw_amp_step{i});
    rp_average{i} = rp_average{i}/length(peak.rp_raw_amp_step{i});

    
end

R.LFP.wavelet.pseudoFreq = pseudoFreq;
% R.LFP.wavelet.coeffs = coeffs;
R.LFP.wavelet.scales = scales;
R.LFP.wavelet.peak = peak;
R.LFP.wavelet.sw_average = sw_average;
R.LFP.wavelet.rp_average = rp_average;

% R.LFP.LFP_ripple_rms = LFP_ripple_rms;
% R.LFP.rms_window_ms = rms_window_ms;

R.LFP.ripple_event = ripple_event;
% R.LFP = rmfield(R.LFP, 'LFP');

end

