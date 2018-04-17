function writeExtInitVHDF5(FID,pop,external_init_V, modify)
% write initial condition for membrane potential 
%            FID: file id for writing data
%            pop: is the number of population
% external_init_V: is the vector of initial membrane potential for all neurons
%                 of this population

if nargin == 3
    modify = 0;
end

if modify == 0
    h5create(FID,['/config/pops/pop',num2str(pop-1),'/SETINITV/external_init_V'], length(external_init_V));
end

h5write(FID,['/config/pops/pop',num2str(pop-1),'/SETINITV/external_init_V'],external_init_V);

end


    
