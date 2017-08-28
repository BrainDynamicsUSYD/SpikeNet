function writePopSeedHDF5(FID, pop_ind, seed, modify)
% manualy seed RNG seed for this neuron popluation
%       FID: file id for writing data
%   pop_ind: neuron population index
%      seed: positive integer
%
% For example, writePopPara(FID, 1, 1234)

if nargin == 3
    modify = 0;
end

% for C/C++ index convetion
pop_ind = pop_ind-1;

% check input
if mod(seed,1) ~= 0 || seed < 0
    disp('The seed should be a positive integer!\n')
else
    if modify == 0
        h5create(FID,['/config/pops/pop',num2str(pop_ind),'/SEED001/seed'], 1);
    end
    h5write(FID,['/config/pops/pop',num2str(pop_ind),'/SEED001/seed'], seed);
end

end
