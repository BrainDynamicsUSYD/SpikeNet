function [R] = get_aCC_pop(R)

fprintf('\t Getting auto-correlation function distribution for each population...\n');

acc_sample_num = 10^3; %10^4; % number of sampling pairs
CC_kernel_width = 40; % ms, kernel length


% Dumpe fields
N = R.N;
spike_hist = R.reduced.spike_hist;
dt = R.reduced.dt;
Num_pop = R.Num_pop;
num_spikes = R.num_spikes;

% define lag vector
lag_max = 100; % ms
lag_inc = 5; % ms

lag_steps = 0:ceil(lag_inc/dt):ceil(lag_max/dt);
lags = lag_steps * dt;


% Define kernel
kernel_type = 'square';
CC_kernel = spike_train_kernel_YG(CC_kernel_width, dt, kernel_type);

% corrcoef within each population
CC_pop = cell(Num_pop,1);
for pop_ind = 1:Num_pop
    if nnz(num_spikes{pop_ind}) > 0
    
        if acc_sample_num >= N(pop_ind)
            samples = 1:N(pop_ind);
        else
            samples = randperm(N(pop_ind), acc_sample_num);
        end
        
        % Calculate pair-wise corrcoef
        for i = 1:length(samples)
            for j = 1:length(lags)
                CC_pop{pop_ind}(i,j) = CorrCoefYG(spike_hist{pop_ind}(samples(i),1:end-lag_steps(j)), ...
                                                  spike_hist{pop_ind}(samples(i),1+lag_steps(j):end), CC_kernel);
            end
            
            % display progress
            if i > 1
                fprintf(repmat('\b',1,10));
            end
            if i < acc_sample_num
                fprintf('%10g',i);
            end
        end
        % Record results
        R.Analysis.aCC_pop_kernel_type = kernel_type;
        R.Analysis.aCC_pop_kernel_width = CC_kernel_width;
        R.Analysis.aCC_pop = CC_pop;
        R.Analysis.aCC_lag = lags;
    end
end

end


