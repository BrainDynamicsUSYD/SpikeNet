function [ R ] = get_SWR( R )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

LFP = R.LFP{1};
dt = R.dt;
fs = 1/(dt*1e-3); % sampling frequency (Hz)

% Butterworth filter
order = 4; % 4th order
lowFreq = 10; % broad band (10?250 Hz)
hiFreq = 250;
Wn = [lowFreq hiFreq]/(fs/2);
[b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.

LFP_broad = filter(b,a,LFP);

% Butterworth filter

order = 4; % 4th order
lowFreq = 80; % ripple band (80?180 Hz)
hiFreq = 180;
Wn = [lowFreq hiFreq]/(fs/2);
[b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.

LFP_ripple = filter(b,a,LFP_broad);

R.SWR.LFP_broad = LFP_broad;
R.SWR.LFP_ripple = LFP_ripple;
end

