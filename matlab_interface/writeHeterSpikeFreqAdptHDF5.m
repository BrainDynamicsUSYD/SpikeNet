function writeHeterSpikeFreqAdptHDF5(FID, pop_ind, dg_K_heter, start_step, end_step, dg_K)
% write heterogeneous spike-frequency adaptation (SFA) setting for specific
% duration: start_step:end_step
% ref: Alessandro Treves, 1993, Mean-field analysis of neuronal spike dynamics
%          FID: file id for writing data
%      pop_ind:
%   dg_K_heter: heterogeneous dg_K for all neurons (uS=miuSiemens)
%start_step, end_step: the heterogeneous SFA works from start step to end step.
%         dg_K: the normal SFA's delta g_k (default value: 0.01 uS)

% for C/C++ index convetion
pop_ind = pop_ind-1;
start_step = start_step - 1;
end_step = end_step - 1;

hdf5write(FID,['/config/pops/pop',num2str(pop_ind),...
    '/INIT010/spike_freq_adpt'],1,'WriteMode','append');

hdf5write(FID,['/config/pops/pop',num2str(pop_ind),...
    '/INIT010/dg_K_heter'],dg_K_heter,'WriteMode','append');

hdf5write(FID,['/config/pops/pop',num2str(pop_ind),...
    '/INIT010/start_step'],start_step,'WriteMode','append');

hdf5write(FID,['/config/pops/pop',num2str(pop_ind),...
    '/INIT010/end_step'],end_step,'WriteMode','append');

if nargin == 6
    hdf5write(FID,['/config/pops/pop',num2str(pop_ind),...
        '/INIT010/dg_K'],dg_K,'WriteMode','append');
end
end