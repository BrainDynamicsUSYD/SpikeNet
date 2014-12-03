function [ R ] = switch_dyanmics_order_para( R )

% order parameter is defined as the infomation rate estimated by shannon
% block entropy of order 2 (H2) of the switch sequence detected.
% 
% order_para = 
%    length(switch sequence) * H2(switch sequence) / simulation_time
%
% Thus the unit for the order parameter is bit per second

% dump fields
Mnum = R.ExplVar.Mnum;
dt = R.dt;
step_tot = R.step_tot;
seq_cell = R.cluster.switch_seq; % a cell

% calculate
T = dt*step_tot/1000; % sec
order_para = zeros(1,length(seq_cell));
for i = 1:length(seq_cell)
    [ H2, ~ ] = symbolic_block_entropy(seq_cell{i}, Mnum); % block entropy of order 2
    num_switches = length(seq_cell{i})-1;
    order_para(i) = H2*num_switches/T; % unit: bit per sec
end


% output result
R.cluster.order_para = order_para;

end