function plot_neuron_sample( R, pop, sample_ind, seg_input )
% this function completely emulates the single neuron behavior in the C++
% simulator


% dump parameters
dt =  R.dt;
step_tot = R.step_tot;
dt_reduced =  R.reduced.dt;
step_tot_reduced = R.reduced.step_tot;

check_cpp_matlab_match = 1;

if nargin == 2
    for sample_ind = 1:length(R.neuron_sample.neuron_ind{pop})
        
        V = R.neuron_sample.V{pop}(sample_ind,:);
        % I_leak = R.neuron_sample.I_leak{pop}(sample_ind,:);
        I_AMPA = R.neuron_sample.I_AMPA{pop}(sample_ind,:);
        I_GABA = R.neuron_sample.I_GABA{pop}(sample_ind,:);
        I_ext = R.neuron_sample.I_ext{pop}(sample_ind,:);
        if isfield( R.neuron_sample, 'I_K')
            I_K = R.neuron_sample.I_K{pop}(sample_ind,:);
        else
            I_K = zeros(size(I_ext));
        end
        neuron_ind = R.neuron_sample.neuron_ind{pop}(sample_ind);
        fprintf('Current sample neuron index = %d.\n', neuron_ind);
        
        spikes = find( R.spike_hist{pop}(neuron_ind,:) );
        
        
        Cm = R.PopPara{pop}.Cm;
        tau_ref = R.PopPara{pop}.tau_ref;
        V_rt = R.PopPara{pop}.V_rt;
        V_lk = R.PopPara{pop}.V_lk;
        g_leak = R.PopPara{pop}.g_lk;
        V_th = R.PopPara{pop}.V_th;
        if isfield(R.PopPara{pop}, 'V_K');
            V_K = R.PopPara{pop}.V_K;
            dg_K = R.PopPara{pop}.dg_K;
            tau_K = R.PopPara{pop}.tau_K;
        else
            V_K = -85.0; % mV
            dg_K = 0.01; %(uS=miuSiemens)
            tau_K = 80; % ms
        end
        
         % simulate the stuff again in matlab (xx_new)
        % reproduce the behavior of the C++ code
        if check_cpp_matlab_match == 1
            spikes_new = [];
            V_new = zeros(1,step_tot);
            V_new(1) = V(1);
            tau_steps = round(tau_ref/dt);
            ref_tmp = 0;
            
            g_K = zeros(1,step_tot);
            I_K = zeros(1,step_tot);
            for t = 1:(step_tot-1)
                % update_spikes
                if ref_tmp == 0 && (V_new(t) >= V_th || (t == 1 && spikes(1) == 1))
                    ref_tmp = tau_steps;
                    spikes_new = [spikes_new t];
                    g_K(t) = g_K(t) + dg_K;
                end
                if any( strcmp(fieldnames(R.ExplVar), 'SpikeFreqAapt') ) && R.ExplVar.SpikeFreqAapt == 1
                    g_K(t) = g_K(t)*exp(-dt/tau_K);
                    I_K(t) = -g_K(t)*(V(t)-V_K);
                    g_K(t+1) = g_K(t);
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
                    
                    V_dot =  (I_leak + I_AMPA(t) + I_GABA(t) + I_ext(t) + + I_K(t))/Cm;
                    
                    %V_dot = (I_leak + I_th)/Cm;
                    
                    V_new(t+1) = V_new(t) + V_dot*dt;
                end
            end
            
            spikes = spikes(:)';
            spikes_new = spikes_new(:)';
            if length(spikes_new) ~= length(spikes) || nnz(sort(spikes_new) - sort(spikes)) > 0
                warning('Simulator does not behave identically as c++ code!');
            end
            
        end
        
        
        % find threshold current for neuron to fire
        I_th = g_leak*(V_th - V_lk); % important!!!!!!!!!!!!!
        I_tot_mean = mean(I_AMPA + I_ext + I_GABA + I_K);
        I_tot_std = std(I_AMPA + I_ext + I_GABA + I_K);
        I_ext_mean = mean(I_ext);
        I_ext_std = std(I_ext);
        I_E_mean = mean(I_AMPA + I_ext);
        I_E_std = std(I_AMPA + I_ext);
        I_I_mean = mean(I_GABA + I_K);
        I_I_std = std(I_GABA + I_K);
        IE_ratio = I_I_mean/I_E_mean;
        
        fprintf('\t I-E ratio = %.3g.\n', IE_ratio);
        fprintf('\t I_E = %.3g (%.3g) including I_ext = %.3g (%.3g).\n', I_E_mean, I_E_std, I_ext_mean, I_ext_std);
        fprintf('\t I_I is %.3g (%.3g).\n',  I_I_mean, I_I_std);
        fprintf('\t I_th = %.3g\n', I_th);
        fprintf('\t I_tot = %.3g (%.3g).\n',  I_tot_mean, I_tot_std);
        
        
        
       
        % Segmetation
        seg_size = 4*10^5; % 40 sec
        seg_size_reduced = round(seg_size*dt/dt_reduced); % be careful!!!
        
        if nargin < 4
            seg_input = 1:ceil(step_tot/seg_size);
        end
        
        for seg = seg_input
            seg_ind = get_seg(step_tot, seg_size, seg);
            seg_ind_reduced = get_seg(step_tot_reduced, seg_size_reduced, seg);
            
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
            
            
            h_ccs = figure('numbertitle','off','name','check_current_stats','color','w', 'position', [680   498   562   600]);
            window_ms = 100; %ms
            window = round(window_ms/dt);
            x = seg_ind*dt*10^-3;
            
            I_E = I_AMPA + I_ext;
            I_I = I_GABA + I_K;
            
            I_E = I_E(seg_ind);
            I_I = I_I(seg_ind);
            
            I_E_std = movingstd(I_E, window);
            I_E_mean = movingmean(I_E, window);
            
            I_I_std = movingstd(I_I, window);
            I_I_mean = movingmean(I_I, window);
            
            I_tot_std = movingstd(I_E+I_I, window);
            
            
            %%%%%%%%%%%%%
            ax(1) = subplot(5,1,1);
            hold on;
            shadedErrorBar(x, I_E_mean, I_E_std, 'r');
            shadedErrorBar(x, I_I_mean, I_I_std, 'b');
            shadedErrorBar(x, I_E_mean+I_I_mean, I_tot_std, 'k');
            plot(x, I_th*ones(size(x)), 'k--');
            ylabel('E/I current mean (std)')
            xlim([min(x) max(x)]);
            
            %%%%%%%%%%%%%
            ax(2) = subplot(5,1,2);
            line([spikes; spikes]*dt*10^-3, [zeros(size(spikes)); ones(size(spikes))]);
            
            num_spikes = nnz(spikes>=min(seg_ind) & spikes<=max(seg_ind));
            ylabel('spikes')
            xlim([min(x) max(x)]);
            
            %%%%%%%%%%%%%
            ax(3) = subplot(5,1,3);
            hold on;
            V_step = 10;
            V_seg = V(seg_ind);
            
            plot(x(1:V_step:end), V_seg(1:V_step:end))
            ymin = min(V_seg);
            ymax = max(V_seg);
            yrange = ymax - ymin;
            ylim([ymin-0.2*yrange,  ymax+0.2*yrange]);
            
            ylabel('membrane potential')
            xlim([min(x) max(x)]);
            
            %%%%%%%%%%%%%
            % moving autocorrelation of V/I
            ax(4) = subplot(5,1,4);
            window_length =  0.5*10^4; % 0.5 sec
            lagNum = 0.1*10^4; % up to 100ms
            sliceNum = 10^3; % the more the better
            [mid_points, lags, acc_mat] = moving_autocorr(I_E, window_length, lagNum, sliceNum);
            imagesc( x(mid_points), lags*dt, acc_mat);
            ylabel('lags (ms)');
            xlim([min(x) max(x)]);
            
            
            %%% individual neuron rate
            ax(5) = subplot(5,1,5);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %     pop_ind = 1;
            %     sigma_gaussian = 50; % ms, which is width???
            %
            %     % Dump fields
            %     N = R.N;
            %
            %     sample_size = 2000; % sample neurons
            %     % down-sampling
            %     if N(pop_ind) >= sample_size
            %         ind_sample = ceil(linspace(1,N(pop_ind),sample_size));
            %     else
            %         ind_sample = 1:1:N(pop_ind);
            %     end
            %
            %     % Dump fields
            %     spike_hist = R.reduced.spike_hist{pop_ind}( ind_sample,  seg_ind_reduced);
            %
            %     % Gaussian filter
            %     kernel = spike_train_kernel_YG(sigma_gaussian, dt_reduced, 'gaussian'); % be careful about the dt here! it's in ms
            %
            %     neuron_rate = zeros(size(spike_hist));
            %     disp('This may take a while...');
            %     for i = 1:length(ind_sample)
            %         neuron_rate(i,:) = SpikeTrainConvolve(spike_hist(i,:), kernel);
            %     end
            %
            %     % imagesc( seg_ind(1:10:end)*dt, 1:length(ind_sample), neuron_rate(:, 1:10:end) );
            %     imagesc( seg_ind_reduced(1:10:end)*dt_reduced*10^-3, 1:length(ind_sample), neuron_rate(:, 1:10:end) );
            
            %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %     % Dump fields
            %     Rate = R.cluster.rate(:,seg_ind_reduced);
            %     Rate = Rate(:,1:10:end); % reduce resolution
            %     cluster_membership = R.cluster.label(neuron_ind);
            %     Mnum = R.ExplVar.Mnum;
            %     imagesc(seg_ind_reduced(1:10:end)*dt_reduced*10^-3, 1:Mnum, Rate);
            %     xlim([min(x) max(x)]);
            
            
            if any( strcmp(fieldnames(R), 'cluster') )
                this_cluster = R.cluster.label(neuron_ind);
                num_spikes_cluster = full(sum( R.reduced.spike_hist{pop}(R.cluster.label == this_cluster, seg_ind_reduced), 1));
                % Plot number of refractory neurons
                x_reduced = seg_ind_reduced*dt_reduced*10^-3;
                line([x_reduced; x_reduced], [zeros(1, length(x_reduced)); num_spikes_cluster], 'color','b');
                
                xlim([min(x) max(x)]);
                
                xlabel(sprintf('t (sec), neuron from cluster %d', this_cluster));
                
                linkaxes(ax,'x');
            end
            
            
            
            %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% The following does not seem to be
            %     % right
            %     figure(3);
            %
            %     freq_range = 1:0.2:100; % Hz
            %     fs = 1/dt*1000; % sampling frequency in Hz
            %     [~,ff,tt,pp] = spectrogram(I_E, window_length,round(window_length*0.9), freq_range, fs, 'yaxis');  %
            %     % Setting 'yaxis' to display frequency on the y-axis and time on the x-axis
            %     % pp is a matrix representing the Power Spectral Density (PSD) of each segment
            %
            %     tt_jump = 1;
            %     tt = tt(1:tt_jump:end);
            %     pp_db = transpose(10*log10(abs(pp(:,1:tt_jump:end))));
            %
            %     subplot(2,1,1)
            %     waterfall(ff, tt, pp_db);
            %     xlabel('Hz');
            %     ylabel('Time (sec)');
            %     zlabel('PSD (dB)');
            %     set(gca,'xscale','log');
            %
            %     subplot(2,1,2)
            %     imagesc(ff, tt, pp_db);
            %     xlabel('Hz');
            %     ylabel('Time (sec)');
            
            
            %
            % lagNum = 10*10^4;
            % x = I_I_mean + I_E_mean;
            % y = I_tot_std;
            %
            % [acfx,lags] = autocorr(x, lagNum);
            % fx = fit(lags', abs(acfx'), 'exp1'); % be very careful about the curve fitting result!
            % figure(4);hold on;
            % plot(lags, acfx);
            % plot(lags, fx.a*exp(fx.b*lags),'r');
            % tau_x = -1/fx.b % the value changes a lot over time!!!!
            %
            %
            % [acfy,lags] = autocorr(y, lagNum);
            % fy = fit(lags', abs(acfy'), 'exp1');
            % figure(5);hold on;
            % plot(lags, acfy);
            % plot(lags, fy.a*exp(fy.b*lags),'r');
            % tau_y = -1/fy.b
            %
            %
            % [xcf,lags] = crosscorr(x,y, lagNum);
            % figure(6);hold on;
            % plot(lags, xcf);
            % fxy = fit(lags((end-1)/2 : end)', -abs(xcf((end-1)/2 : end )'), 'exp1'); % be very careful here!! here I used a dodgy hack
            % plot(lags((end-1)/2 : end), fxy.a*exp(fxy.b*lags((end-1)/2 : end)),'r');
            %
            
            
            next = input(sprintf('\t Next time-frame?'));
            close all;
        end
        next = input('Next sample neuron?');
        close all;
    end
    
    
    
end