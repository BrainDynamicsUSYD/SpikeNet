function Data = ReduceSolution(Data)
    % Re-format the raw spike history data into natural structure with
    % reduced resolution
    % Dump relevant fields
    Num_pop = Data.Num_pop;
    step_tot = Data.step_tot;
    N = Data.N;
    dt = Data.dt;
    spike_hist_compressed = Data.spike_hist_compressed;
    num_spikes = Data.num_spikes;
    num_ref = Data.num_ref;
    
    reduced_dt = 1; % (ms)
    reduced_step_length = round(reduced_dt/dt);
    reduced_step_tot = ceil(step_tot/reduced_step_length);
    reduced_spike_hist = cell(Num_pop,1);
    reduced_num_spikes = cell(Num_pop,1);
    reduced_num_ref = cell(Num_pop,1);
    
    for pop_ind = 1:Num_pop
        %if nnz(num_spikes{pop_ind}) > 0
            reduced_num_spikes_temp = [num_spikes{pop_ind} zeros(1, reduced_step_tot*reduced_step_length-step_tot)]; % padding
            reduced_num_spikes_temp = sum(reshape(reduced_num_spikes_temp, reduced_step_length, reduced_step_tot),1);
            reduced_T_ind_full = cell2mat(arrayfun(@(x, y) repmat(x, [1 y]), 1:reduced_step_tot, reduced_num_spikes_temp, 'UniformOutput', false));
            reduced_spike_hist{pop_ind} = sparse(spike_hist_compressed{pop_ind}, reduced_T_ind_full, true(size(reduced_T_ind_full)), N(pop_ind), reduced_step_tot);
            reduced_num_spikes{pop_ind} = full(sum(reduced_spike_hist{pop_ind},1));
            reduced_num_ref{pop_ind} = num_ref{pop_ind}(round(linspace(1,step_tot,reduced_step_tot))); % down-sampling
        %end
    end

    Data.reduced.dt = reduced_dt;
    Data.reduced.step_tot = reduced_step_tot;
    Data.reduced.spike_hist = reduced_spike_hist;
    Data.reduced.num_spikes = reduced_num_spikes;
    Data.reduced.num_ref = reduced_num_ref;

end