function writeSynFilename(FID, syn_filename)

fprintf(FID, '%s\n', '> SYNF001');
fprintf(FID, '%s\n', syn_filename); % no comma!
fprintf(FID, '\n');

end