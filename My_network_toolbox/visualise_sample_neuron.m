function visualise_sample_neuron(N, dt, V_sample, I_GJ_sample, I_chemical_sample, I_ext_sample, num_spikes_cell)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[Num_pop,step_tot] = size(V_sample);
h_matrix = zeros(5,Num_pop); % matrix for subplot handle

T = (1:step_tot)*dt;

figure('NumberTitle','Off','Name','Activities of sample neurons','position', [50, 300, 900, 668]);

% Manually plot spikes
th = -55; % threshold for firing, make sure it's correct!
spike_peak = 0;
for pop_ind = 1:Num_pop
    for i = 2:step_tot
        if V_sample(pop_ind,i) >= th && V_sample(pop_ind,i-1) < th
            V_sample(pop_ind,i) = spike_peak;
        end
    end
end


for pop_ind = 1:Num_pop
    % Gobal number of spikes
    h_matrix(1,pop_ind) = subplot(5,Num_pop,pop_ind);hold on;
    line([T; T], [zeros(1, step_tot); num_spikes_cell{pop_ind}/N(pop_ind)], 'Color', [255 30 30]/255);
    if pop_ind == 1
        ylabel('% Firing');
    else
        % set(gca,'ytick',[]);
    end
    set(gca,'xtick',[]);
    
    
    % V
    h_matrix(2,pop_ind) = subplot(5,Num_pop,pop_ind+Num_pop*1);hold on;
    plot(T,V_sample(pop_ind,:),'r');
    if pop_ind == 1
        ylabel('Potential (mV)');
    else
        % set(gca,'ytick',[]);
    end
    set(gca,'xtick',[]);
    
    % Chemical current
    h_matrix(3,pop_ind) = subplot(5,Num_pop,pop_ind+Num_pop*2);hold on;
    plot(T,I_chemical_sample(pop_ind,:),'g');
    if pop_ind == 1
        ylabel('I_c_h_e_m_i_c_a_l (?)');
    else
        % set(gca,'ytick',[]);
    end
    set(gca,'xtick',[]);

    
    % External current
    h_matrix(4,pop_ind) = subplot(5,Num_pop,pop_ind+Num_pop*3);hold on;
    plot(T,I_ext_sample(pop_ind,:),'w');
    if pop_ind == 1
        ylabel('I_e_x_t (?)');
    else
        % set(gca,'ytick',[]);
    end
    set(gca,'xtick',[]);

    
    % Electrical current
    h_matrix(5,pop_ind) = subplot(5,Num_pop,pop_ind+Num_pop*4);hold on;
    plot(T,I_GJ_sample(pop_ind,:),'y');
    if pop_ind == 1
        ylabel('I_G_J (?)');
    else
        % set(gca,'ytick',[]);
    end
    xlabel(['t (ms), sampled from pop No.', num2str(pop_ind)]);
    xlim([0, step_tot*dt]); % only one is enough
    
    
end

% Background color
set(h_matrix(:)', 'Color', [0, 0, 0]);

% Link axes to synchronise them when zooming
for pop_ind = 1:Num_pop
    linkaxes(h_matrix(:,pop_ind),'x');
end

% Set initial y-axes to be same, since only one property can be linked
for i = 1:5
    ylimData = [inf -inf]; % initialise ylim
    % find common ylim
    for pop_ind = 1:Num_pop
        ylimData_new = ylim(h_matrix(i,pop_ind));
        ylimData(1) = min(ylimData(1),ylimData_new(1));
        ylimData(2) = max(ylimData(2),ylimData_new(2));
    end
    % set common ylim
    for pop_ind = 1:Num_pop
        ylim(h_matrix(i,pop_ind),ylimData);
    end
end

end








