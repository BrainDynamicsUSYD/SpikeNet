function [R] = get_ISI_low_high(R)
% this is a bloody complicated piece of code!
fprintf('\t Getting separate ISI distributions and CV^2 of ISI for high and low firing states separately...\n');

minimum_spikes = 5; % 3? 5? for CV2_ISI estimation

pop_ind = 1;

% dump fields
N = R.N;
dt = R.reduced.dt;
spike_hist = R.reduced.spike_hist;
label = R.cluster.label;
num_spikes = R.reduced.num_spikes;


seq_1st = R.cluster.sym_seq_1st_theta(2,:); % 2 for 6Hz threshold

c_num = length(unique(label)); % number of clusters

c_high_du = cell(1,c_num);
c_high_start = cell(1, c_num);
c_low_du = cell(1, c_num);
c_low_start = cell(1, c_num);

for i = 1:c_num
    c_high_low = seq_1st;
    c_high_low( c_high_low ~= i ) = 0; % high states of other clusters are still low states for the current cluster
    cut_off = 0; % do not cut head and tail of c_high_low
    [~, high_du_tmp, low_du_tmp, high_start_tmp, low_start_tmp] = seq_postprocess(c_high_low, 1, cut_off); % use dt=1 to get time step
    
    c_high_du{i} = high_du_tmp;
    c_high_start{i} = high_start_tmp;
    c_low_du{i} = low_du_tmp;
    c_low_start{i} = low_start_tmp;
end

%

ISI_cluster_high = []; % ISI 
ISI_cluster_high_ind = []; % neuron index
ISI_cluster_low = [];
ISI_cluster_low_ind = [];


if nnz(num_spikes{pop_ind}) > 0
    for c = 1:c_num

        % high firing state
        for h = 1:length(c_high_du{c})
            a = c_high_start{c}(h);
            b = a + c_high_du{c}(h) - 1;
            for i = find(label == c)
                spike_temp = find(spike_hist{pop_ind}(i, a:b)); % in time step in lieu of ms!!!
                if length(spike_temp) >= 2
                    ISI_tmp = (spike_temp(2:end)-spike_temp(1:end-1))*dt; % in ms
                    ISI_cluster_high = [ISI_cluster_high ISI_tmp];
                    ISI_cluster_high_ind = [ISI_cluster_high_ind i*ones(size(ISI_tmp)) ];
                end
            end
        end
        
        % low firing state
        for h = 1:length(c_low_du{c})
            a = c_low_start{c}(h);
            b = a + c_low_du{c}(h) - 1;
            for i = find(label == c)
                spike_temp = find(spike_hist{pop_ind}(i, a:b)); % in time step in lieu of ms!!!
                if length(spike_temp) >= 2
                    ISI_tmp = (spike_temp(2:end)-spike_temp(1:end-1))*dt; % in ms
                    ISI_cluster_low = [ISI_cluster_low ISI_tmp];
                    ISI_cluster_low_ind = [ISI_cluster_low_ind i*ones(size(ISI_tmp)) ];
                end
            end
        end
        
    end
end


% get CV^2 of ISI for high and low states separately


CV2_ISI_high = zeros(1, N(1) );
CV2_ISI_low = zeros(1, N(1) );

for i = 1:N(1)
    % high
    ISI_tmp = ISI_cluster_high( ISI_cluster_high_ind == i );
    if length(ISI_tmp) >= minimum_spikes
        CV2_ISI_high(i) = (std(ISI_tmp)/mean(ISI_tmp))^2;% CV^2
    else
        CV2_ISI_high(i) = NaN;
    end
    
    % low
    ISI_tmp = ISI_cluster_low( ISI_cluster_low_ind == i );
    if length(ISI_tmp) >= minimum_spikes
        CV2_ISI_low(i) = (std(ISI_tmp)/mean(ISI_tmp))^2;% CV^2
    else
        CV2_ISI_low(i) = NaN;
    end
    
end


% record results
R.Analysis.ISI_cluster_high = ISI_cluster_high;
R.Analysis.ISI_cluster_high_ind = ISI_cluster_high_ind;
R.Analysis.ISI_cluster_low = ISI_cluster_low;
R.Analysis.ISI_cluster_low_ind = ISI_cluster_low_ind;
R.Analysis.CV2_ISI_cluster_high = CV2_ISI_high;
R.Analysis.CV2_ISI_cluster_low = CV2_ISI_low;


end