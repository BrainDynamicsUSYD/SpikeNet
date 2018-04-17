function writeLFPRecordHDF5(FID, pop_ind, LFP_neurons)
%  write data recording for local field potential
%         FID: file id for writing data
%     pop_ind: 
% LFP_neurons: logical vector that specifies which neurons contribute the
%              LFP measure. For multiple LFP measures, use multiple rows.


% for C/C++ index convetion
pop_ind = pop_ind-1;

% write
[n_LFP, N] = size(LFP_neurons);
if n_LFP > N % if the size is not right
    warning('For multiple LFP measures, use multiple rows.')
else
      LFP_neurons = LFP_neurons'; % columns
%     LFP_neurons = LFP_neurons(:); % vectorization
    hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SAMP005/LFP_neurons'],LFP_neurons,'WriteMode','append');
end

end