function [R] = get_rate_spon_evoked(R)
fprintf('\t Getting firing rate of spontanous and evoked activiy, respectively...\n');
% dump fields
dt = R.dt;
step_tot = R.step_tot;
N = R.N;
Num_pop = R.Num_pop;
num_spikes = R.num_spikes;
spike_hist = R.spike_hist;
try
    input_file = h5read([R.stamp(1:end-18),'.h5'],'/config/pops/pop0/file_current_input/fname');
    start_step = h5read(char(input_file),'/start_step');
    end_step = h5read(char(input_file),'/end_step');
    if step_tot <= end_step
        end_step = step_tot;
    end
    step_evoked = end_step - start_step;
    %
    
    rate_spon = cell(Num_pop,1); % Initialise average firing rate of spontaneous activity
    rate_evoked = cell(Num_pop,1); % Initialise average firing rate of evoked activity
    for pop = 1:Num_pop
        if nnz(num_spikes{pop}) > 0
            rate_evoked{pop} = full(sum(spike_hist{pop}(:,start_step:end_step-1),2))/(dt/1000*step_evoked);
            rate_spon{pop} = full(sum(spike_hist{pop}(:,1:start_step),2))/(dt/1000*start_step);% the last part is not pure spontaneous
        end
    end
    % Record results
    R.Analysis.rate_spon = rate_spon;
    R.Analysis.rate_evoked = rate_evoked;
catch
    warning('There is not external current file!');
    % try to get start and end step from explanatory variables
    try
        start_step = get_varible_from_in_h5(R,'sti_start_step');
        end_step = get_varible_from_in_h5(R,'sti_end_step');
        if step_tot <= end_step
            end_step = step_tot;
        end
        step_evoked = end_step - start_step;
        %
        rate_spon = cell(Num_pop,1); % Initialise average firing rate of spontaneous activity
        rate_evoked = cell(Num_pop,1); % Initialise average firing rate of evoked activity
        for pop = 1:Num_pop
            if nnz(num_spikes{pop}) > 0
                rate_evoked{pop} = full(sum(spike_hist{pop}(:,start_step:end_step-1),2))/(dt/1000*step_evoked);
                rate_spon{pop} = full(sum(spike_hist{pop}(:,[1:start_step, end_step:step_tot]),2))/(dt/1000*(step_tot-step_evoked));
            end
        end
        % Record results
        R.Analysis.rate_spon = rate_spon;
        R.Analysis.rate_evoked = rate_evoked;
    catch
        warning('There is not external stimulus!');
    end
end
end
