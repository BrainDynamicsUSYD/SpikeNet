function [popA, indA, popB, indB] = pairs_sample_from_network(N, sample_num)
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

[pairs] = rand_unique_pairs(N_tot, sample_num);

popA = ind_tot2pop(pairs(1,:)); % population index of neuron A
indA = ind_tot2ind(pairs(1,:)); % neuron index of neuron A
popB = ind_tot2pop(pairs(2,:)); % population index of neuron B
indB = ind_tot2ind(pairs(2,:)); % neuron index of neuron B

end