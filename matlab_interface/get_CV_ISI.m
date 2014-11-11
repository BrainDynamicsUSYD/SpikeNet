function [R] = get_CV_ISI(R)

fprintf('\t Getting ISI coefficient of variation...\n');
minimum_spikes = 5; % 3? 5?
% dump fields
dt = R.dt;
N = R.N;
Num_pop = R.Num_pop;
num_spikes = R.num_spikes;
spike_hist = R.spike_hist;
rate = R.Analysis.rate;
ins = 0; % insufficient data
%
CV_ISI = cell(Num_pop,1);% CV_ISI = STD[ISI]/MEAN[ISI], coefficient of variation (note that STD = sqrt(VAR) )
for pop_ind = 1:Num_pop
    if nnz(num_spikes{pop_ind}) > 0
        CV_ISI{pop_ind} = zeros(N(pop_ind),1);
        % warning for insufficient data
        if min(rate{pop_ind}) < minimum_spikes;
            warning('insufficient data for estimating CV_ISI for each neuron!');
            ins = 1;
        end
        % get CV_ISI for neurons with enough spikes
        for i = 1:N(pop_ind)
            spike_temp = (find(spike_hist{pop_ind}(i,:)))*dt; % in (ms) !
            if length(spike_temp) >= minimum_spikes
                ISI_temp = spike_temp(2:end)-spike_temp(1:end-1);
                CV_ISI{pop_ind}(i) = std(ISI_temp)/mean(ISI_temp);% CV
            end
        end
        
    end
end

% CV_ISI averaged over the network
if ins == 0
    CV_ISI_overall = 0;
    for pop_ind = 1:Num_pop
        CV_ISI_overall = CV_ISI_overall+sum(CV_ISI{pop_ind})/sum(N);
    end
    % Record results
    R.Analysis.CV_ISI_overall = CV_ISI_overall;
end


% Record results
R.Analysis.CV_ISI = CV_ISI;


end