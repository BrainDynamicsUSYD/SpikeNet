
function writeInitCond(FID, r_V0, p_fire)
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
    fprintf(FID, '%s\n', '> INIT011');
    fprintf(FID, '%.6f,', r_V0); fprintf(FID,'\n');
    fprintf(FID, '%.6f,', p_fire); fprintf(FID,'\n\n');
end
end

