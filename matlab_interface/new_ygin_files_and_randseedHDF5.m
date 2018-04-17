function [name] = new_ygin_files_and_randseedHDF5(loop_num)


% Creat ygin file
% using loop_num in filename to ensure unique naming!
% Otherwise overwriting may occur when using PBS.
date_now = datestr(now,'yyyymmddHHMM-SSFFF');
name = [ sprintf('%04g-', loop_num), date_now ];
name=[name,'_in.h5'];

fprintf('Data file name is: \n%s\n', name ); % write the file name to stdout and use "grep ygin" to extract it

% seed the matlab rand function! The seed is global.
% Be very careful about that!!!!!!!!!!!
scan_temp = textscan(date_now,'%s','Delimiter','-');
rand_seed = loop_num*10^5+eval(scan_temp{1}{2}); % scan_temp{1}{2} is the SSFFF part
alg='twister';
rng(rand_seed,alg); %  Note that this effect is global!
%NOTE: MATLAB has new h5create and h5write functions but they are less capable (e.g. can't
%write strings, need to create the dataset before writing) so for now use
%the more capable hdf5write function though it may be removed in future
%versions
hdf5write(name,'Original_filename',name);
hdf5write(name,'/MATLAB/rng_seed',rand_seed,'WriteMode','append');
hdf5write(name,'/MATLAB/rng_alg',alg,'WriteMode','append');
% h5create(name,'/MATLAB/rand_alg',1,alg);

fprintf('Random number generator seed is: %f\n', rand_seed );

end