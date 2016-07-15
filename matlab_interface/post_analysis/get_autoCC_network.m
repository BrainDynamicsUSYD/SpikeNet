function [R] = get_autoCC_network(R)
fprintf('\t Getting network activity auto-correlation...\n');

Num_pop = R.Num_pop;
step_tot = R.reduced.step_tot;
dt = R.reduced.dt;
num_spikes = R.reduced.num_spikes;

% network activity auto correlation
num_spikes_tot = zeros(1, step_tot);
for i = 1:Num_pop
    num_spikes_tot = num_spikes_tot + num_spikes{i};
end

max_lag = 50; % ms

[ac,lags] = autocorr( num_spikes_tot, round(max_lag/dt) );

R.Analysis.autoCC_network = ac;
R.Analysis.autoCC_lag = lags*dt;


end

