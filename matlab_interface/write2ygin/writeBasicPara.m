function writeBasicPara(FID, dt, step_tot, N)
% write basic parameters
%      FID: file id for writing data
%       dt: time step size in ms
% step_tot: total number of simulation steps
%        N: vector for number of neurons in each population

%fprintf(FID, '%s\n', '# number of neurons in each population // N1, N2, ...,');
fprintf(FID, '%s\n', '> INIT001');
fprintf(FID, '%d,', N); fprintf(FID,'\n\n');

%fprintf(FID, '%s\n', '# time step length and total number of steps // dt(ms), step_tot,');
fprintf(FID, '%s\n', '> INIT002');
fprintf(FID, '%.9f, %d,\n\n', [dt, step_tot]);


end

