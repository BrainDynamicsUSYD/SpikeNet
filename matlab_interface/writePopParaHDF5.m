function writePopParaHDF5(FID, pop_ind, varargin)
% write non-default neuron model parameters
%       FID: file id for writing data
%   pop_ind: neuron population index
%  varargin: var_name1, var_value1, var_name2, var_value2, ...
%
% For example, writePopPara(FID, 1, "Cm", 1.0, "tau_ref", 3.0)
%
% defualt parameters (using coherent units: msec+mV+nF+miuS+nA)
% Cm = 0.25; % nF
% tau_ref = 2.0; % absolute refractory time (ms), 3.0
% % Potential constants (mV)
% V_rt = -60.0;   % Reset, usually same as leak reversal, use -75.0 to model relative refractory period??
% V_lk = -70.0;   % Leak reversal, -70.0
% V_th = -50.0;   % Threshold // -55.0
% % Leak conductance
% g_lk = 0.0167;  % (miuS), time constants=Cm/gL=15 ms!
% % Reversal potential for external currents
% V_ext = 0.0;      % mV


% for C/C++ index convetion
pop_ind = pop_ind-1;

% check input
var_num = length(varargin)/2;
if mod(var_num,1) ~= 0
    disp('wrong PopPara input format!\n')
    disp('var_num: ');disp(var_num);
else

    para_str = [];
    for i = 1:var_num
        if strcmp(varargin{i*2-1}, 'V_ext')
            warning('Unless you know what you are doing, V_ext should usually be identical to V_ex in SynPara!')
        end
        para_str = [para_str, varargin{i*2-1}, ',' , num2str(varargin{i*2}), ',']; %#ok<AGROW>
    end
    hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/PARA001/para_str_ascii'], double(para_str),'WriteMode','append');
end

end
