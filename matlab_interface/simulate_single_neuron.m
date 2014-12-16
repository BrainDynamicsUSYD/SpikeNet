function [spikes, spikes_new] = simulate_single_neuron( R, pop, sample_ind, seg )
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
        V_new(t+1) = V_new(t) + V_dot*dt;
    end
end

spikes = spikes(:)';
spikes_new = spikes_new(:)';




figure(1);
set(gcf,'color','w');
V_step = 1;

% Segmetation
seg_size = 4*10^4; 
seg_ind = get_seg(step_tot, seg_size, seg);

ax1 = subplot(3,1,1);
line([spikes; spikes], [zeros(size(spikes)); ones(size(spikes))]);
xlim([min(seg_ind) max(seg_ind)]);

ax2 = subplot(3,1,2);
line([spikes_new ; spikes_new ], [zeros(size(spikes_new)); ones(size(spikes_new))]);

x = seg_ind;
V_seg = V(seg_ind);
V_new_seg = V_new(seg_ind);

ax3 = subplot(3,1,3);
hold on;
plot(x(1:V_step:end), V_seg(1:V_step:end), 'b', x(1:V_step:end), V_new_seg(1:V_step:end), 'r')
ymin = min(min(V_seg), min(V_new_seg));
ymax = max(max(V_seg), max(V_new_seg));
yrange = ymax - ymin;
ylim([ymin-0.2*yrange,  ymax+0.2*yrange]);

linkaxes([ax1, ax2, ax3],'x');


end