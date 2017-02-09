function writePopSeedHDF5(FID, pop_ind, seed)
% manualy seed RNG seed for this neuron popluation
%       FID: file id for writing data
%   pop_ind: neuron population index
%      seed: positive integer
%
% For example, writePopPara(FID, 1, 1234)


% for C/C++ index convetion
pop_ind = pop_ind-1;

% check input
if mod(seed,1) ~= 0 || seed < 0
    disp('The seed should be a positive integer!\n')
else
    hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SEED001/seed'], seed,'WriteMode','append');
end

end
