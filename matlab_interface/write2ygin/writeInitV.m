function writeInitV(FID, p_fire)
%fprintf(FID, '%s\n', '# Random initial distributions for membrane potentials // p_fire_pop_1, p_fire_pop_2, ...');
fprintf(FID, '%s\n', '> INIT003');
fprintf(FID, '%.6f,', p_fire);
fprintf(FID,'\n\n');
end

