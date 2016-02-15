function [ R ] = get_stPR( R )
%population coupling measured by the size of spike-triggered population
%rate at zero time lag
%   Calculate population coupling as defined in
%   [1] Michael Okun et al, 2015, Diverse coupling of neurons to populations in sensory cortex

fprintf('\t Getting stPR (may take several minute)...\n');

% spikes detected with 1ms resolution
pop = 1;
spike_hist = R.reduced.spike_hist;
dt = R.reduced.dt; 
N = R.N;

% Gaussian smoothing kernel of half-width (std?) 12 ms
kernel_width = 12; %ms
kernel_type = 'Gaussian_Hz';
sm_kernel = spike_train_kernel_YG(kernel_width, dt, kernel_type);

% population coupling
f_pop = conv( full(sum(double(spike_hist{pop}))), sm_kernel,'same');

max_lag = 500; % 500 ms
lagNum = max_lag/dt;
stPR = zeros(N(1),lagNum*2+1);


for i = 1:N(1)
    f_i = conv( full(double(spike_hist{pop}(i,:))), sm_kernel,'same'); % Hz
    f_i_norm = sum(full(double(spike_hist{pop}(i,:)))); % number of spikes fired (L1-norm?!)
    f_pop_i = f_pop-mean(f_pop) - (f_i-mean(f_i)); % Hz
    [stPR_tmp, lags] = xcorr(f_i, f_pop_i, lagNum);   % eq.A: xcorr(X,Y) --> sum(X.*Y(lag))
    stPR_tmp = stPR_tmp * (dt/1000); % eq.B: convert back to Hz
    stPR(i,:) = stPR_tmp/f_i_norm; % eq.C: normalize
    % the above three equations eq.A-C give the stPR as in the appendix of
    % ref [1]
end

R.stPR.stPR_full = stPR; %in Hz
R.stPR.lags = lags*dt/1000; % sec
R.stPR.kernel_width = kernel_width; %ms
R.stPR.kernel_type = kernel_type;
R.stPR.c = stPR(:,lagNum+1); % size at zero time lag

end

