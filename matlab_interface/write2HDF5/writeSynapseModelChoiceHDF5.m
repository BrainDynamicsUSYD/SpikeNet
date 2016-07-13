function writeSynapseModelChoiceHDF5(FID, model_choice)
%  writeSynapseModelChoice(FID, model_choice)
%   model_choice = 1: See Gu, Yifan, Gong, Pulin, 2016, The dynamics of
%   memory retrieval in hierarchical networks: a modeling study 
%   model_choice = 2: See Keane, A., Gong, P., 2015, Propagating Waves Can
%   Explain Irregular Neural Dynamics 


if model_choice <= 0 || mod(model_choice, 1) ~= 0
    warning('model_choise must be an positive integer!')
elseif model_choice ~= 1 % 1 is default, no need to change it
    model_choice = model_choice-1; % for C/C++ index convetion
    hdf5write(FID,'/config/syns/INIT013/model_choice',model_choice,'WriteMode','append');
end

end