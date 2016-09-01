function writeInitCondHDF5(FID, r_V0, p_fire)
% write initial condition for membrane potential distribution and firing
%    FID: file id for writing data
%   r_V0: set init V distribution to be [V_rt, V_rt + (V_th-V_rt)*r_V0] 
% p_fire: init firing probability
%
% Note that both r_V0 and p_fire should be a vector (defined for all the
% populations).

if max(r_V0) > 1 || max(p_fire) > 1 || min(r_V0) < 0 || min(p_fire) < 0
    error('r_V0 and p_fire must be within 0 and 1!')
else
    for i=1:length(r_V0)
        hdf5write(FID,['/config/pops/pop',num2str(i-1),'/INIT011/r_V0'],r_V0(i),'WriteMode','append');
        hdf5write(FID,['/config/pops/pop',num2str(i-1),'/INIT011/p_fire'],p_fire(i),'WriteMode','append');
    end
end
end

