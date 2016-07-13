function [FID, FID_syn] = new_ygin_files_and_randseed(loop_num, syn)

if nargin == 1
    syn = 1;
end

% Creat ygin file
% using loop_num in filename to ensure unique naming!
% Otherwise overwriting may occur when using PBS.
date_now = datestr(now,'yyyymmddHHMM-SSFFF');
name = [ sprintf('%04g-', loop_num), date_now ];
FID = fopen([name,'.ygin'], 'w'); % creat file
if  syn == 1
    FID_syn = fopen([name,'.ygin_syn'], 'w'); % creat file
else
    FID_syn = -1;
end
fprintf('Data file name is: \n%s\n', strcat(name,'.ygin') ); % write the file name to stdout and use "grep ygin" to extract it

% seed the matlab rand function! The seed is global.
% Be very careful about that!!!!!!!!!!!
scan_temp = textscan(date_now,'%s','Delimiter','-');
rand_seed = loop_num*10^5+eval(scan_temp{1}{2}); % scan_temp{1}{2} is the SSFFF part
rng(rand_seed,'twister'); %  Note that this effect is global!
fprintf('Random number generator seed is: %f\n', rand_seed );

end