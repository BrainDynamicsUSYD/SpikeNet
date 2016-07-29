function [ R ] = get_SWR( R )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% ref: Hippocampal place-cell sequences depict future paths to remembered goals.


LFP = R.LFP.LFP{1};
[no, steps] = size(LFP);
dt = R.dt;
fs = 1/(dt*1e-3); % sampling frequency (Hz)


% Butterworth filter
order = 4; % 4th order
lowFreq = 10; % broad band (10?250 Hz)
hiFreq = 250;
Wn = [lowFreq hiFreq]/(fs/2);
[b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.

for i = 1:no
LFP_broad(i,:) = filter(b,a,LFP(i,:)); %#ok<AGROW>
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
ripple_event.no_std = 1.5;% 3 std above the mean
% Discard transient data
transient_ms = 200; %ms;
ripple_event.transient_ms = transient_ms;
transient_steps = round(transient_ms/dt);
ripple_event.index1 = cell(1,no);
ripple_event.index2 = cell(1,no);
ripple_event.Hz = zeros(1,no);
for i = 1:no
    hil_mean_tmp = mean(LFP_ripple_hilbert(i,transient_steps:end)); % Discard transient data
    hil_std_tmp = std( LFP_ripple_hilbert(i,transient_steps:end)); % Discard transient data
    ripple_seq_tmp = (LFP_ripple_hilbert(i,:) - hil_mean_tmp) > ripple_event.no_std*hil_std_tmp;
    ripple_seq_tmp(1:transient_steps) = 0; % Discard transient data
    ripple_seq_tmp = persistence_requirement( ripple_seq_tmp, ripple_min_steps );
    cutoff = 1;
    [~, ripple_du, flat_du, ripple_start, ~] = seq_postprocess(ripple_seq_tmp, 1, cutoff);
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
    % index3 = zeros(1,N_events);
    for r = 1:N_events
        spike_tmp = spike_hist(LFP_nuerons,ripple_start(r):ripple_start(r)+ripple_du(r)-1);
        ripple_event.index1{i}(r) =  sum(sum(spike_tmp))/spike_tot; % percentage of spikes in each SWR
        participated = (sum(spike_tmp,2))>0;
        ripple_event.index2{i}(participated) = ripple_event.index2{i}(participated) + 1/N_events; % number of neurons participated in each SWR
    end 
end


% Output results
R.LFP.LFP_broad = LFP_broad;
R.LFP.LFP_ripple = LFP_ripple;
R.LFP.LFP_ripple_hilbert = LFP_ripple_hilbert;
R.LFP.lowFreq = lowFreq;
R.LFP.hiFreq = hiFreq;
R.LFP.gauss_width = width;
R.LFP.ripple_event = ripple_event;


end

