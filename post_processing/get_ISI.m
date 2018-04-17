function [R] = get_ISI(R)

minimum_spikes = 5; % 3? 5? for CV^2 calculation
% dump fields
dt = R.dt;
N = R.N;
Num_pop = R.Num_pop;
num_spikes = R.num_spikes;
spike_hist = R.spike_hist;

%
fprintf('\t Getting ISI distribution and CV^2 of ISI...\n');
ISI_lumped = cell(Num_pop,1); % ISI
ISI_lumped_ind = cell(Num_pop,1); % neuron index
CV2_ISI = cell(Num_pop,1);% CV2_ISI = (STD[ISI]/MEAN[ISI])^2, squared coefficient of variation (note that STD = sqrt(VAR) )

for pop_ind = 1:Num_pop
    if nnz(num_spikes{pop_ind}) > 0
        
        % warning for insufficient data for CV^2 of ISI
        CV2_ISI{pop_ind} = zeros(N(pop_ind),1);
        spike_tot = sum(num_spikes{pop_ind},2);
        if min(spike_tot) < minimum_spikes;
            warning('insufficient data for estimating CV_ISI for each neuron!');
        end
        
        
        for i = 1:N(pop_ind)
            spike_temp = find(spike_hist{pop_ind}(i,:)); % in time step in lieu of ms!!!
            % get ISI
            if length(spike_temp) >= 2
                ISI_temp = (spike_temp(2:end)-spike_temp(1:end-1))*dt; % in ms
                ISI_lumped{pop_ind} = [ISI_lumped{pop_ind} ISI_temp];
                ISI_lumped_ind{pop_ind} = [ISI_lumped_ind{pop_ind} i*ones(size(ISI_temp)) ];
            end
            % get CV^2 of ISI
            if length(spike_temp) >= minimum_spikes
                ISI_temp = (spike_temp(2:end)-spike_temp(1:end-1))*dt; % in ms
                CV2_ISI{pop_ind}(i) = (std(ISI_temp)/mean(ISI_temp))^2;% CV^2
            else
                CV2_ISI{pop_ind}(i) = NaN;
            end
        end
    end
end

% record results
R.Analysis.ISI_lumped = ISI_lumped;
R.Analysis.ISI_lumped_ind = ISI_lumped_ind;
R.Analysis.CV2_ISI = CV2_ISI;



end