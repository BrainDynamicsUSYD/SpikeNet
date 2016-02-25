
function writeInitV(FID, p_fire)
% write initial condition for membrane potential distribution
%    FID: file id for writing data
% p_fire: initial firing rate for each population (vector, value 0~1)
%
% Note that the membrane potential distributions are uniform and decided by
% the given p_fire

%fprintf(FID, '%s\n', '# Random initial distributions for membrane potentials // p_fire_pop_1, p_fire_pop_2, ...');
fprintf(FID, '%s\n', '> INIT003');
fprintf(FID, '%.6f,', p_fire);
fprintf(FID,'\n\n');
end

