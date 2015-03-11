function [R] = get_CC_network(R)
fprintf('\t Getting correlation coefficient distribution (network-wide sampling)...\n');
corrcoef_sample_num = 10^3; %10^4; % number of sampling pairs
CC_kernel_width = 40; % ms, kernel length

% Dumpe fields
N = R.N;
spike_hist = R.reduced.spike_hist;
dt = R.reduced.dt;
Num_pop = R.Num_pop;
num_spikes = R.num_spikes;

% Define kernel
kernel_type = 'square';
CC_kernel = spike_train_kernel_YG(CC_kernel_width, dt, kernel_type);

% sample pairs from entire network (all the populations)
NNZ_product = 1;
for pop_ind = 1:Num_pop
    NNZ_product = NNZ_product * nnz(num_spikes{pop_ind}); % if any one of them is zero
end

if NNZ_product > 0
    [popA, indA, popB, indB] = pairs_sample_from_network(N,corrcoef_sample_num);
    corrcoef_sample_num = length(popA); % Note that length(popA) <= corrcoef_sample_num !!!!!!!!!!!!!
    % Calculate pair-wise corrcoef
    CC_sample = zeros(1,corrcoef_sample_num);
    for i = 1:corrcoef_sample_num
        CC_sample(i) = CorrCoefYG(spike_hist{popA(i)}(indA(i),:), spike_hist{popB(i)}(indB(i),:), CC_kernel);
        % display progress
        if i > 1
            fprintf(repmat('\b',1,10));
        end
        if i < corrcoef_sample_num
            fprintf('%10g',i);
        end
    end
    % Record results
    R.Analysis.CC_network_kernel_type = kernel_type;
    R.Analysis.CC_network_kernel_width = CC_kernel_width;
    R.Analysis.CC_network = CC_sample;
end
end







