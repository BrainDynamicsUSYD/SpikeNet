function [ R ] = get_stPR( R )
%population coupling measured by the size of spike-triggered population
%rate at zero time lag
%   Calculate population coupling as defined in
%   [1] Michael Okun et al, 2015, Diverse coupling of neurons to populations in sensory cortex

fprintf('\t Getting stPR (may take several minute)...\n');
pop = 1;

% sample a sub-set of neurons from the entire population
sample_num = 100;
sample_ind = randperm(R.N(pop), sample_num);

% spikes detected with 1ms resolution
spike_hist = R.reduced.spike_hist{pop}(sample_ind,:);
dt = R.reduced.dt; 

% Gaussian smoothing kernel of half-width (std?) 12 ms
kernel_width = 12; %ms
kernel_type = 'Gaussian_Hz';
sm_kernel = spike_train_kernel_YG(kernel_width, dt, kernel_type);

% population coupling
f_pop = conv( full(sum(double(spike_hist))), sm_kernel,'same');

max_lag = 500; % 500 ms
lagNum = max_lag/dt;
stPR = zeros(sample_num,lagNum*2+1);

for i = 1:sample_num
    f_i = conv( full(double(spike_hist(i,:))), sm_kernel,'same'); % Hz
    f_i_norm = sum(full(double(spike_hist(i,:)))); % number of spikes fired (L1-norm?!)
    f_pop_i = f_pop-mean(f_pop) - (f_i-mean(f_i)); % Hz
    [stPR_tmp, lags] = xcorr(f_i, f_pop_i, lagNum);   % eq.A: xcorr(X,Y) --> sum(X.*Y(lag))
    stPR_tmp = stPR_tmp * (dt/1000); % eq.B: convert back to Hz
    stPR(i,:) = stPR_tmp/f_i_norm; % eq.C: normalize
    % the above three equations eq.A-C give the stPR as in the appendix of
    % ref [1]
end


%%%%%%%%%%%%%% shuffling
% Implementation of the Raster Marginals Model, as introduced in
% "Population rate dynamics and multineuron firing patterns in 
%  sensory cortex", Journal of Neuroscience.
spike_hist_shuffle = spike_hist';
shuffle_num = 30*nchoosek(sample_num,2); 
c = ceil(rand(shuffle_num,2)*sample_num); % two randomly selected columns
for i = 1:shuffle_num
  I = spike_hist_shuffle(:,c(i,1)) + spike_hist_shuffle(:,c(i,2)) == 1; % where the 2 columns don't coincide
  cA = spike_hist_shuffle(I,[c(i,1) c(i,2)]); % a copy of the part that matters, to make it run faster
  i01 = find(cA(:,1) == 0);
  i10 = find(cA(:,1) == 1);
  toFlip = ceil(min(length(i01), length(i10))/2); % how many 01s & 10s to flip
  i01 = i01(randperm(length(i01)));
  i01 = i01(1:toFlip);
  i10 = i10(randperm(length(i10)));
  i10 = i10(1:toFlip);
  % the flip itself:
  cA(i01,1) = true; cA(i01,2) = false;
  cA(i10,1) = false; cA(i10,2) = true;
  spike_hist_shuffle(I,[c(i,1) c(i,2)]) = cA;  
end;
spike_hist_shuffle = spike_hist_shuffle'; % the shuffled spike history
% calculate stPR for shuffled spikes
f_pop = conv( full(sum(double(spike_hist_shuffle))), sm_kernel,'same');
stPR_shuffle = zeros(sample_num,lagNum*2+1);
for i = 1:sample_num
    f_i = conv( full(double(spike_hist_shuffle(i,:))), sm_kernel,'same'); % Hz
    f_i_norm = sum(full(double(spike_hist_shuffle(i,:)))); % number of spikes fired (L1-norm?!)
    f_pop_i = f_pop-mean(f_pop) - (f_i-mean(f_i)); % Hz
    [stPR_tmp, lags] = xcorr(f_i, f_pop_i, lagNum);   % eq.A: xcorr(X,Y) --> sum(X.*Y(lag))
    stPR_tmp = stPR_tmp * (dt/1000); % eq.B: convert back to Hz
    stPR_shuffle(i,:) = stPR_tmp/f_i_norm; % eq.C: normalize
end

R.stPR.sample_ind = sample_ind;
R.stPR.stPR_full = stPR; %in Hz
R.stPR.lags = lags*dt/1000; % sec
R.stPR.kernel_width = kernel_width; %ms
R.stPR.kernel_type = kernel_type;
R.stPR.c = stPR(:,lagNum+1); % size at zero time lag

R.stPR.stPR_full_shuffle =  stPR_shuffle;
R.stPR.c_shuffle =  stPR_shuffle(:,lagNum+1);
R.stPR.c_norm = R.stPR.c/nanmedian(R.stPR.c_shuffle); % normlize by the median value from shuffled data

end

