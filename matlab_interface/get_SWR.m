function [ R ] = get_SWR( R )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% ref: Hippocampal place-cell sequences depict future paths to remembered goals.


LFP = R.LFP.LFP{1};
[no, steps] = size(LFP);
dt = R.dt;
fs = 1/(dt*1e-3); % sampling frequency (Hz)
width = 12.5; %ms
[ Kernel ] = spike_train_kernel_YG( width, dt, 'gaussian_unit' );

%
transient_ms = 0; %ms;
transient_steps = round(transient_ms/dt);
if steps <= transient_steps
    error('Not enough LFP data!')
end

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
    LFP_ripple_hilbert(i,:) = conv(abs(hilbert(LFP_ripple(i,:))), Kernel,'same'); %#ok<AGROW>
end


% Output results
R.LFP.transient_steps = transient_steps;
R.LFP.LFP{1} = R.LFP.LFP{1}(:,transient_steps+1:end);
R.LFP.LFP_broad = LFP_broad(:,transient_steps+1:end);
R.LFP.LFP_ripple = LFP_ripple(:,transient_steps+1:end);
R.LFP.LFP_ripple_hilbert = LFP_ripple_hilbert(:,transient_steps+1:end);
R.LFP.lowFreq = lowFreq;
R.LFP.hiFreq = hiFreq;
R.LFP.gauss_width = width;
end

