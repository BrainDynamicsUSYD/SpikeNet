function [ R ] = get_stPR( R, varargin )
%population coupling measured by the size of spike-triggered population
%rate at zero time lag
%   Calculate population coupling as defined in
%   [1] Michael Okun et al, 2015, Diverse coupling of neurons to populations in sensory cortex

fprintf('\t Getting stPR (may take several minute)...\n');
pop = 1;

n_trials = 10;
sample_num = 66;
n_sample_region = 2000;

for i = 1:length(varargin)/2
    var_name = varargin{2*i-1};
    var_value = varargin{2*i};
     if isnumeric(var_value)
        eval([var_name, '=', num2str(var_value), ';']);
     else
         eval([var_name, '=''', var_value, ''';']);
     end
end

dt = R.reduced.dt;
max_lag = 500; % 500 ms
lagNum = max_lag/dt;

% Gaussian smoothing kernel of "half-width" 12/sqrt(2) ms
% FWHM = sigma * sqrt(8*log(2))
FWHM = 12/sqrt(2)*2;
kernel_std = FWHM / sqrt(8*log(2)) ; %ms
kernel_type = 'Gaussian_Hz';
sm_kernel = spike_train_kernel_YG(kernel_std, dt, kernel_type);

% sample a sub-set of neurons from the entire population
fw = sqrt(R.N(pop));
hw = (fw-1)/2;
if mod(hw,1) ~= 0
    error('Not a supported grid.');
end

rad = sqrt(n_sample_region/pi); % 25 neurons is about 125 miu-meter
Lat = lattice_nD(2,hw);
ind_within = find(Lat(:,1).^2+ Lat(:,2).^2 <= rad^2);


sample_ind = [];
stPR = [];
stPR_shuffle = [];
c_norm = [];
c_shuffle_norm = [];

for n = 1:n_trials
    sample_ind_tmp = ind_within( randperm(length(ind_within), sample_num) );
    
    sample_ind = [sample_ind; sample_ind_tmp(:)]; %#ok<*AGROW>
    % spikes detected with 1ms resolution
    spike_hist = R.reduced.spike_hist{pop}(sample_ind_tmp,:);
   
    
    % population coupling
    f_pop = conv( full(sum(double(spike_hist))), sm_kernel,'same');
    
    stPR_tmp = zeros(sample_num,lagNum*2+1);
    
    for i = 1:sample_num
        f_i = conv( full(double(spike_hist(i,:))), sm_kernel,'same'); % Hz
        f_i_norm = sum(full(double(spike_hist(i,:)))); % number of spikes fired (L1-norm?!)
        f_pop_i = f_pop-mean(f_pop) - (f_i-mean(f_i)); % Hz
        [stPR_tmp_i, lags] = xcorr(f_i, f_pop_i, lagNum);   % eq.A: xcorr(X,Y) --> sum(X.*Y(lag))
        stPR_tmp_i = stPR_tmp_i * (dt/1000); % eq.B: convert back to Hz
        stPR_tmp(i,:) = stPR_tmp_i/f_i_norm; % eq.C: normalize
        % the above three equations eq.A-C give the stPR as in the appendix of
        % ref [1]
    end
    stPR = [stPR; stPR_tmp];
    
%     %%%%%%%%%%%%%% shuffling
%     % Implementation of the Raster Marginals Model, as introduced in
%     % "Population rate dynamics and multineuron firing patterns in
%     %  sensory cortex", Journal of Neuroscience.
%     spike_hist_shuffle = spike_hist';
%     shuffle_num = 100*nchoosek(sample_num,2);
%     c = ceil(rand(shuffle_num,2)*sample_num); % two randomly selected columns
%     for i = 1:shuffle_num
%         I = spike_hist_shuffle(:,c(i,1)) + spike_hist_shuffle(:,c(i,2)) == 1; % where the 2 columns don't coincide
%         cA = spike_hist_shuffle(I,[c(i,1) c(i,2)]); % a copy of the part that matters, to make it run faster
%         i01 = find(cA(:,1) == 0);
%         i10 = find(cA(:,1) == 1);
%         toFlip = ceil(min(length(i01), length(i10))/2); % how many 01s & 10s to flip
%         i01 = i01(randperm(length(i01)));
%         i01 = i01(1:toFlip);
%         i10 = i10(randperm(length(i10)));
%         i10 = i10(1:toFlip);
%         % the flip itself:
%         cA(i01,1) = true; cA(i01,2) = false;
%         cA(i10,1) = false; cA(i10,2) = true;
%         spike_hist_shuffle(I,[c(i,1) c(i,2)]) = cA;
%     end;
%     spike_hist_shuffle = spike_hist_shuffle'; % the shuffled spike history

    spike_hist_shuffle = raster_marginals_shuffling(spike_hist);

    % calculate stPR for shuffled spikes
    f_pop = conv( full(sum(double(spike_hist_shuffle))), sm_kernel,'same');
    stPR_shuffle_tmp = zeros(sample_num,lagNum*2+1);
    for i = 1:sample_num
        f_i = conv( full(double(spike_hist_shuffle(i,:))), sm_kernel,'same'); % Hz
        f_i_norm = sum(full(double(spike_hist_shuffle(i,:)))); % number of spikes fired (L1-norm?!)
        f_pop_i = f_pop-mean(f_pop) - (f_i-mean(f_i)); % Hz
        [stPR_s_tmp, lags] = xcorr(f_i, f_pop_i, lagNum);   % eq.A: xcorr(X,Y) --> sum(X.*Y(lag))
        stPR_s_tmp = stPR_s_tmp * (dt/1000); % eq.B: convert back to Hz
        stPR_shuffle_tmp(i,:) = stPR_s_tmp/f_i_norm; % eq.C: normalize
    end
    
    stPR_shuffle = [stPR_shuffle; stPR_shuffle_tmp];
    
    c_norm_tmp = stPR_tmp(:,lagNum+1)/nanmedian(stPR_shuffle_tmp(:,lagNum+1));
    c_shuffle_norm_tmp = stPR_shuffle_tmp(:,lagNum+1)/nanmedian(stPR_shuffle_tmp(:,lagNum+1));
    
    c_norm = [c_norm; c_norm_tmp(:)];
    c_shuffle_norm = [c_shuffle_norm; c_shuffle_norm_tmp(:)];
end


R.stPR.n_trials = n_trials;
R.stPR.sample_ind = sample_ind;

R.stPR.lags = lags*dt/1000; % sec
R.stPR.kernel_width = kernel_std; %ms
R.stPR.kernel_type = kernel_type;

R.stPR.stPR_full = stPR'; %in Hz
R.stPR.stPR_full_shuffle =  stPR_shuffle';

R.stPR.c_shuffle_norm =  c_shuffle_norm;
R.stPR.c_norm = c_norm; % normlize by the median value from shuffled data

end

