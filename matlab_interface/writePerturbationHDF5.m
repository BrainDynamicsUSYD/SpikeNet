function writePerturbationHDF5(FID, pop_ind, step_perturb)
% write perturbation setting
%          FID: file id for writing data
%      pop_ind:
% step_perturb: the step where one spike is removed (if there is no spike at this step, then the next step)


% for C/C++ index convetion
pop_ind = pop_ind-1;

% write
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT007/step_perturb'],step_perturb,'WriteMode','append');

end

