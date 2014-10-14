function Data = ReformatSpikeHistory(Data)
    % Re-format the raw (compressed) spike history data into natural structure
    % (logical sparse matrix)
    % Dump relevant fields
    Num_pop = Data.Num_pop;
    step_tot = Data.step_tot;
    N = Data.N;
    spike_hist = Data.spike_hist;
    num_spikes = Data.num_spikes;
    
    % Re-format the raw spike history data into natural structure
    spike_hist_natural = cell(Num_pop,1);
    for pop_ind = 1:Num_pop
        if nnz(num_spikes{pop_ind}) > 0
            T_ind_full = cell2mat(arrayfun(@(x, y) repmat(x, [1 y]), 1:step_tot, num_spikes{pop_ind}, 'UniformOutput', false));
            spike_hist_natural{pop_ind} = sparse(spike_hist{pop_ind}, T_ind_full, true(size(T_ind_full)), N(pop_ind), step_tot);
        else
            spike_hist_natural{pop_ind} = sparse(false(N(pop_ind), step_tot));
        end
    end
    Data.spike_hist_compressed = Data.spike_hist;
    Data.spike_hist = spike_hist_natural; % overwrite spike history data (logical sparse matrix)
end