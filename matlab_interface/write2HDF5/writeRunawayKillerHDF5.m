function writeRunawayKillerHDF5(FID, pop_ind, min_ms, runaway_Hz, Hz_ms)
% write runaway activity killer
%            FID: file id for writing data
%         min_ms: minimum simulation time before it can be killed
%     runaway_Hz: the threshold defines runaway activity
%          Hz_ms: the time window over which the populational firing rate
%                 is estimated and compared with runaway_Hz
%
% So if any neuron population's activity exceeds the threshold define by
% (runaway_Hz,Hz_ms), the simulation will be killed.
% Note that only those populations with (by defualt) >100 neurons will be
% considered.

% for C/C++ index convetion
pop_ind = pop_ind-1;


%fprintf(FID, '%s\n', '# runaway killer setting //runaway_steps, runaway_mean_num_ref(0~1),');
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/KILL001/min_ms'],min_ms,'WriteMode','append');
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/KILL001/runaway_Hz'],runaway_Hz,'WriteMode','append');
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/KILL001/Hz_ms'],Hz_ms,'WriteMode','append');
end