
function [fname_E,fname_I]= Adam_stimulus_mean_curr(mean_curr,sigma2,centerx,centery)
% Adam's stimulu: a Guassian shape stimulus
% transfer it to external current
% Adam_stimulus(amplitude,10,168,149)

netsz=250; % MAKE THIS GRIDSIZE OF EXCITATORY NEURONS
neurons_E=1:netsz*netsz*0.75;
neurons_I=1:netsz*netsz*0.25;
Grid = zeros(netsz);
Grid(2:2:netsz,2:2:netsz) = 1;
ind_E = Grid == 0;
ind_I = Grid == 1;

FileName = 'Guassian_Adam';

[x,y] = ndgrid(1:netsz);
d2 = (x-centery).^2 + (y-centerx).^2; % meshgrid is better for avoiding reversal x and y
frame = exp(-d2./(2.*sigma2));
frame = frame/mean(frame(:));

frame_E = frame(ind_E);
frame_I = frame(ind_I);
curr_E=reshape(frame_E,1,netsz*netsz*0.75);% turn the 2D image into a 1D list
curr_I=reshape(frame_I,1,netsz*netsz*0.25);
end_step = 0.8*10^5;
start_step = 0.2*10^5;
% mean_curr=1; % change this value to control the mean current input
% write out the data to a file
% Each frame is a row and consists of a list of neuron indexes that will
% recieve an external current from the file, and a list of currents for
% each of these neurons. Need a minimum of two frames, hence the repeated
% [curr;curr] and [neurons;neurons]
fname_E=[FileName,'_mean_curr',num2str(mean_curr),'_',num2str(centerx),'x',num2str(centery),'_E.h5']; % ***** change the filename if using a different dataset
fname_I=[FileName,'_mean_curr',num2str(mean_curr),'_',num2str(centerx),'x',num2str(centery),'_I.h5'];
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

