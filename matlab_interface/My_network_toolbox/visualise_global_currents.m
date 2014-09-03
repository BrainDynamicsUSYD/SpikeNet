function visualise_global_currents( N, dt, I_chemical_stat, I_GJ_stat, I_ext_stat, num_spikes_cell)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


[~,step_tot,Num_pop] = size(I_chemical_stat); % [mean;stddev]
h_matrix = zeros(4,Num_pop); % matrix for subplot handle

figure('NumberTitle','Off','Name','Statistics of driving currents','position', [50, 300, 900, 668]);

T = (1:step_tot)*dt;

for pop_ind = 1:Num_pop
    % Gobal number of spikes
    h_matrix(1,pop_ind) = subplot(4,Num_pop,pop_ind);hold on;
    line([T; T], [zeros(1, step_tot); num_spikes_cell{pop_ind}/N(pop_ind)], 'Color', [255 30 30]/255);
    if pop_ind == 1
        ylabel('% Firing');
    else
        % set(gca,'ytick',[]);
    end
    set(gca,'xtick',[]);
    
    
    % I_chemical
    h_matrix(2,pop_ind) = subplot(4,Num_pop,pop_ind+Num_pop*1);hold on;
    plot(T,I_chemical_stat(1,:,pop_ind),'g',T,I_chemical_stat(1,:,pop_ind)+I_chemical_stat(2,:,pop_ind),'g:',T,I_chemical_stat(1,:,pop_ind)-I_chemical_stat(2,:,pop_ind),'g:');
    if pop_ind == 1
        ylabel('I_c_h_e_m_i_c_a_l (?)');
    else
        % set(gca,'ytick',[]);
    end
    set(gca,'xtick',[]);
    
    % I_ext
    h_matrix(3,pop_ind) = subplot(4,Num_pop,pop_ind+Num_pop*2);hold on;
    plot(T,I_ext_stat(1,:,pop_ind),'w',T,I_ext_stat(1,:,pop_ind)+I_ext_stat(2,:,pop_ind),'w:',T,I_ext_stat(1,:,pop_ind)-I_ext_stat(2,:,pop_ind),'w:');
    if pop_ind == 1
        ylabel('I_e_x_t (?)');
    else
        % set(gca,'ytick',[]);
    end
    set(gca,'xtick',[]);

    
    % I_GJ
    h_matrix(4,pop_ind) = subplot(4,Num_pop,pop_ind+Num_pop*3);hold on;
    plot(T,I_GJ_stat(1,:,pop_ind),'y',T,I_GJ_stat(1,:,pop_ind)+I_GJ_stat(2,:,pop_ind),'y:',T,I_GJ_stat(1,:,pop_ind)-I_GJ_stat(2,:,pop_ind),'y:');
    if pop_ind == 1
        ylabel('I_G_J (?)');
    else
        % set(gca,'ytick',[]);
    end
    xlabel(['t (ms), from pop No.', num2str(pop_ind)]);
    xlim([0, step_tot*dt]); % only one is enough
    
    
end

% Background color
set(h_matrix(:)', 'Color', [0, 0, 0]);

% Link axes to synchronise them when zooming
for pop_ind = 1:Num_pop
    linkaxes(h_matrix(:,pop_ind),'x');
end

% Set initial y-axes to be same, since only one property can be linked
for i = 1:4
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

