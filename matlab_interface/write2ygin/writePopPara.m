function writePopPara(FID, pop_ind, varargin)
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
% g_lk = 0.0167;   % (miuS), time constants=Cm/gL=15 ms!


% for C/C++ index convetion
pop_ind = pop_ind-1;

% check input
var_num = length(varargin)/2;
if mod(var_num,1) ~= 0
    disp('wrong PopPara input format!\n')
    disp('var_num: ');disp(var_num);
else
    % fprintf(FID, '%s\n', '# non-default neuron population parameter setting, data following: pop_ind, number of parameters,; var_name1,value1,;var_name2,value2,;...');
    fprintf(FID, '%s\n', '> PARA001');
    fprintf(FID, '%d,%d,\n', pop_ind, var_num);
    for i = 1:var_num
        fprintf(FID, '%s,%.9f,\n', varargin{i*2-1}, varargin{i*2});
    end
    fprintf(FID,'\n');
end
end

