function writePopCovRecordHDF5(FID, pop_ind, time_start, time_end)

% for C/C++ index convetion
pop_ind = pop_ind-1;
time_start = time_start - 1;
time_end = time_end - 1;

% write
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SAMP103/time_start'],time_start,'WriteMode','append');
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SAMP103/time_end'],time_end,'WriteMode','append');

end