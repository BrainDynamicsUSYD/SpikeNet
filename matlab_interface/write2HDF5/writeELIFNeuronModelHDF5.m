function writeELIFNeuronModelHDF5(FID,pop,ELIF_VT,ELIF_delT)
%       FID: file id for writing data
%       pop: the population number
%       type: type of neuron model (0=LIF, 1 = Exponential LIF)
%       ELIF_VT: Exponential integrate and Fire threshold
%       ELIF_delT: Exponential integrate and FIre spiking slope factor



type=1; %this specifies elif
    
    
hdf5write(FID,['/config/pops/pop',num2str(pop-1),'/neuron_model'],int32(type),'WriteMode','append');

if type==1
    hdf5write(FID,['/config/pops/pop',num2str(pop-1),'/ELIF/ELIF_VT'],ELIF_VT,'WriteMode','append'); 
    hdf5write(FID,['/config/pops/pop',num2str(pop-1),'/ELIF/ELIF_delT'],ELIF_delT,'WriteMode','append');
end
    
