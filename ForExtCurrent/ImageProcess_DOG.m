
function [fname_E,fname_I]= ImageProcess_DOG(mean_curr,gender,num)
% process the nature image
% filtered by DOG(difference of Gaussian) then make them as
% external current and transfer it to HDF5 file for E and I populations
% This is for bubbles paper's database.
% gender = ['w','m'];num = 1:8

path = 'Add your dataset path here';

% pyr = cell(1,16); pind = cell(1,16);
netsz=250; % MAKE THIS GRIDSIZE OF EXCITATORY NEURONS
neurons_E=1:netsz*netsz*0.75;
neurons_I=1:netsz*netsz*0.25;
Grid = zeros(netsz);
Grid(2:2:netsz,2:2:netsz) = 1;
ind_E = Grid == 0;
ind_I = Grid == 1;

FaceName = [gender,num2str(num),'N'];
face = imread([path,FaceName,'.JPG']);
imwrite(face,['original_',FaceName,'.jpg'])

if length(size(face,3)) == 3
    face = rgb2gray(face);
end

face=imresize(face,[netsz netsz]); % resize the square image to the same size as the excitatory neural population
% get local variance
h1 = imgaussfilt(face,1);
h2 = imgaussfilt(face,2);
face_DOG = h1-h2;
figure('NumberTitle', 'off','Name','face_LV', 'visible', 'off');
imshow(face_DOG,[])
% save the real img (without border)
f = getframe;
real_img = rgb2gray(f.cdata);
imwrite(real_img,['DOG_',FaceName,'.jpg'])
frame = double(face_DOG);
frame_E = frame(ind_E);
frame_I = frame(ind_I);
curr = reshape(frame,1,netsz*netsz);
curr_E=reshape(frame_E,1,netsz*netsz*0.75);% turn the 2D image into a 1D list
curr_I=reshape(frame_I,1,netsz*netsz*0.25);
curr_E=curr_E/mean(curr);  % normalize the image/currents to have mean=1
curr_I=curr_I/mean(curr);
end_step = 0.8*10^6;
start_step = 0.2*10^5;

% write out the data to a file
% Each frame is a row and consists of a list of neuron indexes that will
% recieve an external current from the file, and a list of currents for
% each of these neurons. Need a minimum of two frames, hence the repeated
% [curr;curr] and [neurons;neurons]
fname_E=[FaceName,'_',num2str(mean_curr),'_',num2str(netsz),'x',num2str(netsz),'_LV_E.h5']; % ***** change the filename if using a different dataset
fname_I=[FaceName,'_',num2str(mean_curr),'_',num2str(netsz),'x',num2str(netsz),'_LV_I.h5'];
hdf5write(fname_E,'/current',[curr_E;curr_E]');
hdf5write(fname_E,'/neurons',int32([neurons_E;neurons_E]'),'WriteMode','append');
hdf5write(fname_E,'/frame_rate',1,'WriteMode','append');
hdf5write(fname_E,'/mean_curr',mean_curr,'WriteMode','append');
hdf5write(fname_E,'/end_step',end_step,'WriteMode','append');
hdf5write(fname_E,'/start_step',start_step,'WriteMode','append');
hdf5write(fname_I,'/current',[curr_I;curr_I]');
hdf5write(fname_I,'/neurons',int32([neurons_I;neurons_I]'),'WriteMode','append');
hdf5write(fname_I,'/frame_rate',1,'WriteMode','append');
hdf5write(fname_I,'/mean_curr',mean_curr,'WriteMode','append');
hdf5write(fname_I,'/end_step',end_step,'WriteMode','append');
hdf5write(fname_I,'/start_step',start_step,'WriteMode','append');
end

