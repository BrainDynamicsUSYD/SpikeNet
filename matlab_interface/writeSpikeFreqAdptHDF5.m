function writeSpikeFreqAdptHDF5(FID, pop_ind, dg_K)
% write spike-frequency adaptation setting
% ref: Alessandro Treves, 1993, Mean-field analysis of neuronal spike dynamics
%          FID: file id for writing data
%      pop_ind:
%         dg_K: default value 0.01 (uS=miuSiemens)

% for C/C++ index convetion
pop_ind = pop_ind-1;


hdf5write(FID,['/config/pops/pop',num2str(pop_ind),...
    '/INIT010/spike_freq_adpt'],1,'WriteMode','append');

if nargin == 3
    hdf5write(FID,['/config/pops/pop',num2str(pop_ind),...
        '/INIT010/dg_K'],dg_K,'WriteMode','append');
end

end

