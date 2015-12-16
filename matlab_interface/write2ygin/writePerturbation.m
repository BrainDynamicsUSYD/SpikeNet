function writePerturbation(FID, pop_ind, step_perturb)
% write perturbation setting
%          FID: file id for writing data
%      pop_ind:
% step_perturb: the step where one spike is removed (if there is no spike at this step, then the next step)

fprintf(FID, '%s\n', '> INIT007');
fprintf(FID, '%d,%d,', pop_ind-1, step_perturb);
fprintf(FID,'\n\n');
end

