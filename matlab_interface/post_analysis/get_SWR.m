function [ R ] = get_SWR( R )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% ref: Hippocampal place-cell sequences depict future paths to remembered goals.


LFP = R.LFP.LFP{1};
[no, steps] = size(LFP);

no = 1;

dt = R.dt;
fs = 1/(dt*1e-3); % sampling frequency (Hz)


% Butterworth filter
order = 4; % 4th order
lowFreq = 1; % broad band (1-1250 ? Hz)
hiFreq = 1250;
Wn = [lowFreq hiFreq]/(fs/2);
[b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.

for i = 1:no
    LFP_broad(i,:) = filter(b,a,LFP(i,:)); %#ok<AGROW>
end

% Butterworth filter
order = 4; % 4th order
lowFreq = 1; % broad band (10?250 Hz)
hiFreq = 50;
Wn = [lowFreq hiFreq]/(fs/2);
[b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.

for i = 1:no
    LFP_sharpwave(i,:) = filter(b,a,LFP(i,:)); %#ok<AGROW>
end


% Butterworth filter
order = 4; % 4th order
lowFreq = 150; % ripple band (80?180 Hz)
hiFreq = 250;
Wn = [lowFreq hiFreq]/(fs/2);
[b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
for i = 1:no
    LFP_ripple(i,:) = filter(b,a,LFP(i,:)); %#ok<AGROW>
    % hilbert transformation & gaussian smoothing
    width = 5; %ms
    [ Kernel ] = spike_train_kernel_YG( width, dt, 'gaussian_unit' );
    LFP_ripple_hilbert(i,:) = conv(abs(hilbert(LFP_ripple(i,:))), Kernel,'same'); %#ok<AGROW>
    
%     rms_window_ms = 17; %ms
%     window_steps = round(rms_window_ms/dt);
%     Kernel_rms = ones(1,window_steps)/window_steps;
%     LFP_ripple_rms(i,:) = sqrt(conv(LFP_ripple(i,:).^2, Kernel_rms,'same')); %#ok<AGROW>
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
spike_hist = R.spike_hist{1};
spike_tot = sum(sum(spike_hist));
ripple_event.ripple_min_ms = 15;
ripple_min_steps = round(ripple_event.ripple_min_ms/dt);
ripple_event.no_std = 2;% 3 std above the mean
ripple_event.peak_no_std = 4;% 3 std above the mean
% Discard transient data
transient_ms = 200; %ms;
ripple_event.transient_ms = transient_ms;
transient_steps = round(transient_ms/dt);
ripple_event.index1 = cell(1,no);
ripple_event.index2 = cell(1,no);
ripple_event.Hz = zeros(1,no);
ripple_event.is_SWR = zeros(size(LFP_ripple_hilbert));

iter_num = 100;
hil_mean_baseline_hist = zeros(no,iter_num);
hil_std_baseline_hist = zeros(no,iter_num);
for i = 1:no
    LFP_ripple_hilbert_tmp = LFP_ripple_hilbert(i,:); 
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
    % three indexes for population synchrony
    % ref: Preconfigured, skewed distribution of firing rates in the hippocampus and entorhinal cortex
    LFP_nuerons = logical(R.LFP.LFP_neurons{1}(i,:));
    N_events = length(ripple_du);
    ripple_event.index1{i} = zeros(1,N_events);
    ripple_event.index2{i} = zeros(1,sum(LFP_nuerons));
    ripple_event.index3{i} = zeros(1,sum(LFP_nuerons));
    % index3 = zeros(1,N_events);
    for r = 1:N_events
        spike_tmp = spike_hist(LFP_nuerons,ripple_start(r):ripple_start(r)+ripple_du(r)-1);
        ripple_event.index1{i}(r) =  sum(sum(spike_tmp))/spike_tot; % percentage of spikes in each SWR
        participated = (sum(spike_tmp,2))>0;
        ripple_event.index2{i}(participated) = ripple_event.index2{i}(participated) + 1/N_events; % number of neurons participated in each SWR
        ripple_event.index3{i} = ripple_event.index3{i} + (sum(spike_tmp,2)/N_events)'; % mean number of spikes per SWR
    end
end
ripple_event.hil_mean_baseline = hil_mean_baseline_hist;
ripple_event.hil_std_baseline = hil_std_baseline_hist;



%
freqrange = [R.LFP.lowFreq R.LFP.hiFreq];
Fs = 1000/dt;
fc = centfrq('cmor1.5-1');
scalerange = fc./(freqrange*(1/Fs));
scales = scalerange(end):0.5:scalerange(1);
pseudoFreq = scal2frq(scales,'cmor1.5-1',1/Fs); % pseudo-frequencies
coeffs = cell(no,1);
peak = cell(no,1);
for i = 1:no
    x_tmp = R.LFP.LFP_ripple(i,:);
    coeffs_tmp = abs(cwt(x_tmp,scales,'cmor1.5-1'))';
    peak{i}.freq = [];
    peak{i}.step = [];
    peak{i}.coeff = [];
    
    for j =  1:length(ripple_event.ripple_du_steps{i})
        a_tmp = ripple_event.ripple_start_steps{i}(j);
        l_tmp = ripple_event.ripple_du_steps{i}(j);
        
        [max_coffs_tmp, freq_ind_tmp] = max(coeffs_tmp(a_tmp:a_tmp+l_tmp-1, :), [], 2);
        [max_max_coffs_tmp, ind_tmp] = max(max_coffs_tmp);
        peak{i}.freq = [peak{i}.freq pseudoFreq(freq_ind_tmp(ind_tmp))];
        peak{i}.step = [peak{i}.step a_tmp+ind_tmp-1];
        peak{i}.coeff = [peak{i}.coeff max_max_coffs_tmp];
    end
     
    coeffs{no} = coeffs_tmp;
end

R.LFP.wavelet.pseudoFreq = pseudoFreq;
R.LFP.wavelet.coeffs = coeffs;
R.LFP.wavelet.peak = peak;

% Output results
R.LFP.LFP_broad = LFP_broad;
R.LFP.LFP_ripple = LFP_ripple;
R.LFP.LFP_sharpwave = LFP_sharpwave;
R.LFP.LFP_ripple_hilbert = LFP_ripple_hilbert;
% R.LFP.LFP_ripple_rms = LFP_ripple_rms;
% R.LFP.rms_window_ms = rms_window_ms;
R.LFP.lowFreq = lowFreq;
R.LFP.hiFreq = hiFreq;
R.LFP.gauss_width = width;
R.LFP.ripple_event = ripple_event;


end

