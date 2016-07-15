function [R] = get_crosscorr_network(R)
% according to "A Master Equation Formalism for MAcroscopic Modeling of
% Asynchronous Irregular Activity States, P48, P50
% this a network-wide measure of synchrony
% The value at zero time lag can be used as a convenient single value
% measure with 5ms bin and 500 disjoint pairs (Kumar et al 2008)

fprintf('\t Getting average cross-correlation (network-wide sampling)...\n');
sample_num = 500; % number of sampling pairs
CC_kernel_width = 5; % ms, kernel length
kernel_type = 'square';
% max_lag for auto-correlation
max_lag = 20; %ms


% Dumpe fields
N = R.N;
spike_hist = R.reduced.spike_hist;
dt = R.reduced.dt;
Num_pop = R.Num_pop;
num_spikes = R.num_spikes;

% Define kernel
CC_kernel = spike_train_kernel_YG(CC_kernel_width, dt, kernel_type);


% sample pairs from entire network (all the populations)
NNZ_product = 1;
for pop_ind = 1:Num_pop
    NNZ_product = NNZ_product * nnz(num_spikes{pop_ind}); % if any one of them is zero
end

if NNZ_product > 0
    [popA, indA, popB, indB] = pairs_sample_from_network(N,sample_num);
    sample_num = length(popA); % Note that length(popA) <= corrcoef_sample_num !!!!!!!!!!!!!
    % Calculate pair-wise cross-correlation
    CC_sample = [];
    
    for i = 1:sample_num
        train_A = SpikeTrainConvolve(spike_hist{popA(i)}(indA(i),:),CC_kernel);
        train_B = SpikeTrainConvolve(spike_hist{popB(i)}(indB(i),:),CC_kernel);
        if i == 1
        [CC_tmp, lags] = crosscorr( train_A, train_B, round(max_lag/dt) );
        else
             CC_tmp = crosscorr( train_A, train_B, round(max_lag/dt) );
        end
        CC_sample = [CC_sample; CC_tmp];
        
        % display progress
        if i > 1
            fprintf(repmat('\b',1,10));
        end
        if i < sample_num
            fprintf('%10g',i);
        end
    end
    % Record results
    
    CC_full = mean(CC_sample, 1);
    
    R.Analysis.crosscorr_network_kernel_type = kernel_type;
    R.Analysis.crosscorr_network_kernel_width = CC_kernel_width;
    R.Analysis.crosscorr_network_full = CC_full;
    R.Analysis.crosscorr_network_lags = lags*dt;
    R.Analysis.crosscorr_network_0 = CC_full(lags == 0);
end

end

