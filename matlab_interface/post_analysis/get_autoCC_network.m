function [R] = get_autoCC_network(R)
fprintf('\t Getting network activity auto-correlation...\n');

Num_pop = R.Num_pop;
% step_tot = R.reduced.step_tot;
dt = R.reduced.dt;
num_spikes = R.reduced.num_spikes;

max_lag = min(1000, round( (R.step_tot*R.dt)/ 4) ); % ms

% network activity auto correlation
ac = [];
for i = 1:Num_pop
    [ac_tmp,lags] = autocorr( num_spikes{i}, round(max_lag/dt) );
    ac = [ac;ac_tmp(:)']; %#ok<AGROW>
end

R.Analysis.autoCC_network = ac;
R.Analysis.autoCC_lag = lags*dt;
end

