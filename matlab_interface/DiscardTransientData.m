function Data = DiscardTransientData(Data)
    % Discard trainsient data
    if any( strcmp(fieldnames(Data), 'ExplVar') ) && ...
            ~isempty(Data.ExplVar) && ...
            any( strcmp(fieldnames(Data.ExplVar), 'discard_transient') ) && ...
            Data.ExplVar.discard_transient > 0 % input check
        
        % discard_transient = Data.ExplVar.discard_transient;
        % step_tot_dis = round(Data.step_tot*(1-discard_transient));
        
        discard_transient = Data.ExplVar.discard_transient; % in (ms)
        dt = Data.dt;
        step_transient = round( discard_transient/dt ); % in steps
        step_tot_dis = Data.step_tot - step_transient;
        
        if step_tot_dis <= 0
            warning('discard_transient cannot be larger than total simulation time!');
            return;
        end
        
        % spiking history
        for pop_ind = 1:Data.Num_pop
                Data.spike_hist{pop_ind} = Data.spike_hist{pop_ind}(:,end-step_tot_dis+1:end);
                if nnz(Data.num_spikes{pop_ind}) > 0
                    Data.spike_hist_compressed{pop_ind} = Data.spike_hist_compressed{pop_ind}(sum(Data.num_spikes{pop_ind}(1:end-step_tot_dis))+1:end);
                end
                Data.num_spikes{pop_ind} = Data.num_spikes{pop_ind}(end-step_tot_dis+1:end);
                Data.num_ref{pop_ind} = Data.num_ref{pop_ind}(end-step_tot_dis+1:end);
        end
        
        % neuron sample
        if any( strcmp(fieldnames(Data), 'neuron_sample') ) && ~isempty(Data.neuron_sample)
            fname_cell = fieldnames(Data.neuron_sample);
            for pop_ind = 1:length(Data.neuron_sample.neuron_ind)
                if ~isempty(Data.neuron_sample.neuron_ind{pop_ind})
                    t_dis = find(Data.neuron_sample.t_ind{pop_ind} <= step_transient);
                    Data.neuron_sample.t_ind{pop_ind}(t_dis) = [];
                    Data.neuron_sample.t_ind{pop_ind} = Data.neuron_sample.t_ind{pop_ind} - step_transient;
                    for f = 1:length(fname_cell)
                        fn = fname_cell{f};
                        if ~strcmp(fn,'neuron_ind') && ~strcmp(fn,'t_ind') % not "t_ind" or "neuron_ind" 
                            Data.neuron_sample.(fn){pop_ind}(:,t_dis) = [];
                        end
                    end
                end
            end
        end
        
        %  
        if any( strcmp(fieldnames(Data), 'pop_stats') ) && ~isempty(Data.pop_stats)
            t_dis = 1:step_transient;
            for pop = 1:length(Data.pop_stats.V_mean)
                if ~isempty(Data.pop_stats.V_mean{pop})
                    
                    Data.pop_stats.V_mean{pop}(t_dis) = [];
                    Data.pop_stats.V_std{pop}(t_dis) = [];
                    
                    Data.pop_stats.I_input_mean{pop}(t_dis) = [];
                    Data.pop_stats.I_input_std{pop}(t_dis) = [];
                end
            end
        end
        
        if any( strcmp(fieldnames(Data), 'syn_stats') ) && ~isempty(Data.syn_stats)
            t_dis = 1:step_transient;
            for syn = 1:length(Data.syn_stats)
                Data.syn_stats{syn}.I_mean(t_dis) = [];
                Data.syn_stats{syn}.I_std(t_dis) = [];
            end
        end
        
                
        
        % synapse sample
        if any( strcmp(fieldnames(Data), 'syn_sample') ) && ~isempty(Data.syn_sample)
            for syn_statsnd = 1:length(Data.syn_sample)
                t_dis = find(Data.syn_sample{syn_statsnd}.t_ind <= step_transient);
                Data.syn_sample{syn_statsnd}.t_ind(t_dis) = [];
                Data.syn_sample{syn_statsnd}.t_ind = Data.syn_sample{syn_statsnd}.t_ind - step_transient;
                Data.syn_sample{syn_statsnd}.I(:,t_dis) = [];
                
            end
        end

        % update step_tot
        Data.step_tot = step_tot_dis;
        if Data.step_killed > 0
            Data.step_killed = Data.step_killed - step_transient;
        end
        
    end
end