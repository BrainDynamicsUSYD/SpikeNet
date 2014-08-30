function writeSynPara(FID, varargin)
% write non-default synapse parameters
%       FID: file id for writing data
%  varargin: var_name1, var_value1, var_name2, var_value2, ...
%
% For example, writeSynPara(FID, "V_ex", 1.0, "V_in", -85.0)
%
% defualt parameters (using coherent units: msec+mV+nF+miuS+nA)
% % Potential constants (mV)
% V_ex = 0.0;     % Excitatory reversal, 0.0
% V_in = -80.0;   % Inhibitory reversal, -80.0
% % time-evolution of post-synaptic conductance change
% Dt_trans_AMPA = 1.0; % 0.5
% Dt_trans_GABA = 1.0; % 1.0
% Dt_trans_NMDA = 5.0; % 5.0
% tau_decay_AMPA = 5.0; % 3.0
% tau_decay_GABA = 5.0; % 7.0
% tau_decay_NMDA = 80.0; % 80.0


% check input
var_num = length(varargin)/2;
if mod(var_num,1) ~= 0
    disp('wrong SynPara input format!\n')
    disp('var_num: ');disp(var_num);
else
    % fprintf(FID, '%s\n', '# synapse parameter setting, data following: number of parameters,; var_name1,value1,;var_name2,value2,;...');
    fprintf(FID, '%s\n', '> PARA002');
    fprintf(FID, '%d,\n', var_num);
    for i = 1:var_num
        fprintf(FID, '%s,%.9f,\n', varargin{i*2-1}, varargin{i*2});
    end
    fprintf(FID,'\n');
end
end
