function writeInitVHDF5(FID, p_fire)
% write initial condition for membrane potential distribution
%    FID: file id for writing data
% p_fire: initial firing rate for each population (vector, value 0~1)
%
% Note that the membrane potential distributions are uniform and decided by
% the given p_fire

for i=1:length(p_fire)
    hdf5write(FID,['/config/pops/pop',num2str(i-1),'/INIT003/p_fire'],p_fire(i),'WriteMode','append')
end
   
end

