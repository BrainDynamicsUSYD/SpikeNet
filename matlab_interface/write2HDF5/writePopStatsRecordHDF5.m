function writePopStatsRecordHDF5(FID, pop_ind)

% for C/C++ index convetion
pop_ind = pop_ind-1;
% write
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SAMP003/record'],1,'WriteMode','append');

end