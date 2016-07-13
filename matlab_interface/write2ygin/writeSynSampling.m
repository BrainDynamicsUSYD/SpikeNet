function writeSynSampling(FID, pop_ind_pre, pop_ind_post, syn_type, sample_ind, time_index)
% write data sampling settings for individual neurons
%        FID: file id for writing data
% sample_ind: vector of post-synaptic neuron indices to be sampled in pop_ind
% time_index: logical vector to define time points to be sampled
% 
% Note that time_index should have the length of total simulation steps.
% For example, if step_tot = 10 and the last half time points are to be
% sampled, use time_index = [0,0,0,0,0,1,1,1,1,1]
% 4,000 neurons x 10,000 time points will be 300MB data


if sum((time_index ~= 0) & (time_index ~= 1)) ~= 0
    error('sample_time_index must be logical vectors: ');
end

% for C/C++ index convetion
pop_ind_pre = pop_ind_pre-1;
pop_ind_post = pop_ind_post-1;
syn_type = syn_type-1;
sample_ind = sample_ind-1;


% write
% fprintf(FID, '%s\n', '# synapse data sampling');
fprintf(FID, '%s\n', '> SAMP002');
fprintf(FID, '%d, %d, %d,\n', pop_ind_pre, pop_ind_post, syn_type);
fprintf(FID, '%d,', sample_ind); fprintf(FID,'\n');
fprintf(FID, '%d,', time_index); fprintf(FID,'\n\n');


end
