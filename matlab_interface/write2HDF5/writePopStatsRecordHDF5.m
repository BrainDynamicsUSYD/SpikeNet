function writePopStatsRecordHDF5(FID, pop_ind)

% for C/C++ index convetion
pop_ind = pop_ind-1;

% write

hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SAMP003'],1,'WriteMode','append');
%+++NOT SURE IF THIS IS THE CORRECT PLACE TO PUT THIS???
end