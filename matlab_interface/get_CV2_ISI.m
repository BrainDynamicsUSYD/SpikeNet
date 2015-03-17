function [R] = get_CV2_ISI(R)

fprintf('\t Getting squared coefficient of variation of ISI...\n');
minimum_spikes = 5; % 3? 5?
% dump fields
dt = R.dt;
N = R.N;
Num_pop = R.Num_pop;
num_spikes = R.num_spikes;
spike_hist = R.spike_hist;

%
CV2_ISI = cell(Num_pop,1);% CV2_ISI = (STD[ISI]/MEAN[ISI])^2, squared coefficient of variation (note that STD = sqrt(VAR) )
for pop_ind = 1:Num_pop
    if nnz(num_spikes{pop_ind}) > 0
        CV2_ISI{pop_ind} = zeros(N(pop_ind),1);
        spike_tot = sum(num_spikes{pop_ind},2);
        % warning for insufficient data
        if min(spike_tot) < minimum_spikes;
            warning('insufficient data for estimating CV_ISI for each neuron!');
        end
        % get CV_ISI for neurons with enough spikes
        for i = 1:N(pop_ind)
            spike_temp = (find(spike_hist{pop_ind}(i,:)))*dt; % in (ms) !
            if length(spike_temp) >= minimum_spikes
                ISI_temp = spike_temp(2:end)-spike_temp(1:end-1);
                CV2_ISI{pop_ind}(i) = (std(ISI_temp)/mean(ISI_temp))^2;% CV^2
            else
                CV2_ISI{pop_ind}(i) = NaN;
            end
        end
        
    end
end

% Record results
R.Analysis.CV2_ISI = CV2_ISI;


end