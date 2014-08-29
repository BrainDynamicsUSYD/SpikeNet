function writeRespVar(FID, var_name)
fprintf(FID, '%s\n', '> response variable');
fprintf(FID, '%s,\n', var_name);
end

