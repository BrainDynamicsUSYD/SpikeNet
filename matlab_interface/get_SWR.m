function [ R ] = get_SWR( R )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% ref: Hippocampal place-cell sequences depict future paths to remembered goals.

LFP = R.LFP{1};
dt = R.dt;
fs = 1/(dt*1e-3); % sampling frequency (Hz)
width = 12.5; %ms
[ Kernel ] = spike_train_kernel_YG( width, dt, 'gaussian_unit' );

% Butterworth filter
order = 4; % 4th order
lowFreq = 10; % broad band (10?250 Hz)
hiFreq = 250;
Wn = [lowFreq hiFreq]/(fs/2);
[b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
LFP_broad = filter(b,a,LFP);


% Butterworth filter
order = 4; % 4th order
lowFreq = 150; % ripple band (80?180 Hz)
hiFreq = 250;
Wn = [lowFreq hiFreq]/(fs/2);
[b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.

LFP_ripple = filter(b,a,LFP);

% hilbert transformation & gaussian smoothing
LFP_ripple_hilbert = conv(abs(hilbert(LFP_ripple)), Kernel,'same');


R.SWR.LFP_broad = LFP_broad;
R.SWR.LFP_ripple = LFP_ripple;
R.SWR.LFP_ripple_hilbert = LFP_ripple_hilbert;
R.SWR.lowFreq = lowFreq;
R.SWR.hiFreq = hiFreq;
R.SWR.gauss_width = width;
end

