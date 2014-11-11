function [R] = get_ISI(R)
% dump fields
dt = R.dt;
N = R.N;
Num_pop = R.Num_pop;
num_spikes = R.num_spikes;
spike_hist = R.spike_hist;
rate = R.Analysis.rate;
%
fprintf('\t Getting ISI distribution...\n');
ISI_dist = cell(Num_pop,1); % ISI histogram
for pop_ind = 1:Num_pop
    if nnz(num_spikes{pop_ind}) > 0
        for i = 1:N(pop_ind)
            if rate{pop_ind}(i) >= 2
                spike_temp = find(spike_hist{pop_ind}(i,:)); % in time step in lieu of ms!!!
                Dt_temp = (spike_temp(2:end)-spike_temp(1:end-1))*dt; % in ms
                ISI_dist{pop_ind} = [ISI_dist{pop_ind} Dt_temp];
            end
        end
    end
end
% record results
R.Analysis.ISI_dist = ISI_dist;

end