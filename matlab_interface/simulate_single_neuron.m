function [spikes, spikes_new] = simulate_single_neuron( R, pop, sample_ind, seg_input )
% this function completely emulates the single neuron behavior in the C++
% simulator


% dump parameters
dt =  R.dt;
step_tot = R.step_tot;

V = R.neuron_sample.V{pop}(sample_ind,:);
% I_leak = R.neuron_sample.I_leak{pop}(sample_ind,:);
I_AMPA = R.neuron_sample.I_AMPA{pop}(sample_ind,:);
I_GABA = R.neuron_sample.I_GABA{pop}(sample_ind,:);
I_ext = R.neuron_sample.I_ext{pop}(sample_ind,:);
neuron_ind = R.neuron_sample.neuron_ind{pop}(sample_ind);

spikes = find( R.spike_hist{pop}(neuron_ind,:) );


Cm = R.PopPara{pop}.Cm;
tau_ref = R.PopPara{pop}.tau_ref;
V_rt = R.PopPara{pop}.V_rt;
V_lk = R.PopPara{pop}.V_lk;
g_leak = R.PopPara{pop}.g_lk;
V_th = R.PopPara{pop}.V_th;

% find threshold current for neuron to fire
I_th = g_leak*(V_th - V_lk) % important!!!!!!!!!!!!!
I_E_mean = mean(I_AMPA + I_ext)
I_I_std = std(I_AMPA + I_ext)
I_I_mean = mean(I_GABA)
I_E_std = std(I_GABA)
I_tot_mean = mean(I_AMPA + I_ext + I_GABA)
I_tot_std = std(I_AMPA + I_ext + I_GABA)
EI_ratio = mean(I_AMPA + I_ext)/mean(I_GABA)


% simulate the stuff again in matlab (xx_new)
% reproduce the behavoir of the C++ code

spikes_new = [];
V_new = zeros(1,step_tot);
V_new(1) = V(1);

tau_steps = round(tau_ref/dt);
ref_tmp = 0;

for t = 1:(step_tot-1)
    % update_spikes
    if ref_tmp == 0 && V_new(t) >= V_th
        ref_tmp = tau_steps;
        spikes_new = [spikes_new t];
    end
    if ref_tmp > 0
        ref_tmp = ref_tmp - 1;
        V_new(t) = V_rt;
    end
    % update_V
    if ref_tmp == 0
        I_leak = -g_leak*(V_new(t) - V_lk);
        
        %I_GABA(t) = 0;
        %I_AMPA(t) = 0;
        %I_ext(t) = 0;
        
        V_dot =  (I_leak + I_AMPA(t) + I_GABA(t) + I_ext(t))/Cm;
        
        %V_dot = (I_leak + I_th)/Cm;
        
        V_new(t+1) = V_new(t) + V_dot*dt;
    end
end





spikes = spikes(:)';
spikes_new = spikes_new(:)';


% Segmetation
seg_size = 4*10^5; 

if nargin < 4
    seg_input = 1:ceil(step_tot/seg_size);
end

for seg = seg_input
    seg_ind = get_seg(step_tot, seg_size, seg);
    
    
    % figure('numbertitle','off','name','check_simulation_correctness','color','w');
    % V_step = 1;
    %
    % ax1 = subplot(3,1,1);
    % line([spikes; spikes], [zeros(size(spikes)); ones(size(spikes))]);
    % xlim([min(seg_ind) max(seg_ind)]);
    %
    % ax2 = subplot(3,1,2);
    % line([spikes_new ; spikes_new ], [zeros(size(spikes_new)); ones(size(spikes_new))]);
    %
    % x = seg_ind;
    % V_seg = V(seg_ind);
    % V_new_seg = V_new(seg_ind);
    %
    % ax3 = subplot(3,1,3);
    % hold on;
    % plot(x(1:V_step:end), V_seg(1:V_step:end), 'b', x(1:V_step:end), V_new_seg(1:V_step:end), 'r')
    % ymin = min(min(V_seg), min(V_new_seg));
    % ymax = max(max(V_seg), max(V_new_seg));
    % yrange = ymax - ymin;
    % ylim([ymin-0.2*yrange,  ymax+0.2*yrange]);
    %
    % linkaxes([ax1, ax2, ax3],'x');
    
    
    h_ccs = figure('numbertitle','off','name','check_current_stats','color','w');
    window_ms = 100; %ms
    window = round(window_ms/dt);
    x = seg_ind;
    
    I_E = I_AMPA+I_ext;
    I_I = I_GABA;
    
    I_E = I_E(seg_ind);
    I_I = I_I(seg_ind);
    
    I_E_std = movingstd(I_E, window);
    I_E_mean = movingmean(I_E, window);
    
    I_I_std = movingstd(I_I, window);
    I_I_mean = movingmean(I_I, window);
    
    ax1 = subplot(4,1,1);
    plot(x, I_E_mean, 'r', x, I_I_mean, 'b')
    xlim([min(seg_ind) max(seg_ind)]);
    ylabel('E/I current mean')
    
    ax2 = subplot(4,1,2);
    plot(x, I_E_std, 'r', x, I_I_std, 'b')
    ylabel('E/I current std')
    
    
    ax3 = subplot(4,1,3);
    line([spikes; spikes], [zeros(size(spikes)); ones(size(spikes))]);
    xlim([min(seg_ind) max(seg_ind)]);
    num_spikes = nnz(spikes>=min(seg_ind) & spikes<=max(seg_ind))
    ylabel('spikes')
    
    ax4 = subplot(4,1,4);
    hold on;
    V_step = 1;
    V_seg = V(seg_ind);
    
    plot(x(1:V_step:end), V_seg(1:V_step:end))
    ymin = min(V_seg);
    ymax = max(V_seg);
    yrange = ymax - ymin;
    ylim([ymin-0.2*yrange,  ymax+0.2*yrange]);
    ylabel('membrane potential')
    
    linkaxes([ax1, ax2, ax3, ax4],'x');
    
    
    next = input('next?');
    close(h_ccs);
end



end