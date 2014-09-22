function [R] = get_CC_network(R)
fprintf('\t Getting correlation coefficient distribution (network-wide sampling)...\n');
corrcoef_sample_num = 10^4; % number of sampling pairs
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



% helper function
function [popA, indA, popB, indB] = pairs_sample_from_network(N, corrcoef_sample_num)
% e.g.
% N = [4000, 1000]; sample_num = 10^4
% popA = [1 2 1 1 1 2    1    2 ....] (10^4 elements)
% indA = [1 1 2 3 4 1000 4000 8 ...] (10^4 elements)

N_tot = sum(N);
ind_tot2pop = []; % total index to pop id 1:Num_pop
ind_tot2ind = []; % total index to pop index 1:N(i)

for i = 1:length(N)
    ind_tot2pop = [ind_tot2pop i*ones(1,N(i))];
    ind_tot2ind = [ind_tot2ind 1:N(i)];
end

[pairs] = rand_unique_pairs(N_tot, corrcoef_sample_num);

popA = ind_tot2pop(pairs(1,:)); % population index of neuron A
indA = ind_tot2ind(pairs(1,:)); % neuron index of neuron A
popB = ind_tot2pop(pairs(2,:)); % population index of neuron B
indB = ind_tot2ind(pairs(2,:)); % neuron index of neuron B

end



