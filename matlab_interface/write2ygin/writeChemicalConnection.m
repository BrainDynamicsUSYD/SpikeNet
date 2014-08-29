function writeChemicalConnection(FID, type, i_pre, j_post, I, J, K, D)
% for C/C++ index convetion
type = type -1;
i_pre = i_pre -1;
j_post = j_post -1;
I = I-1; 
J = J-1;
if nargin == 7
    D = zeros(size(I)); % auto zero delay
end
% write (note: no white space!!!)
%fprintf(FID, '%s\n', '# chemical connection // (type,i_pre,j_post);I;J;K;D;');
fprintf(FID, '%s\n', '> INIT006');
fprintf(FID, '%d,', [type;i_pre;j_post]); fprintf(FID,'\n');
fprintf(FID, '%d,', I); fprintf(FID,'\n');
fprintf(FID, '%d,', J); fprintf(FID,'\n');
fprintf(FID, '%.9f,', K); fprintf(FID,'\n'); % note precision, K in miuSiemens
fprintf(FID, '%.3f,', D); fprintf(FID,'\n'); % note precision, D in msec
fprintf(FID, '\n');

end