function writePopStatsRecordHDF5(FID, pop_ind, varargin)

% for C/C++ index convetion
pop_ind = pop_ind-1;

if isempty(varargin)
    % write
    hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SAMP003/record'],1,'WriteMode','append');
elseif length(varargin) == 2
    time_start = varargin{1} - 1;
    time_end = varargin{2} - 1;
    % write
    hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SAMP103/time_start'],time_start,'WriteMode','append');
    hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SAMP103/time_end'],time_end,'WriteMode','append');
end
end