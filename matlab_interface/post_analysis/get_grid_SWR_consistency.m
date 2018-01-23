function R = get_grid_SWR_consistency(R)

% Load data
dt = R.dt;
stamp = R.stamp;
samp_file = [stamp(1:end-3) '0_neurosamp'];
load(['./', samp_file], 'peak','ripple_fit', 'ripple_fit_goodness', 'peak_LFP','RPI', 'ripple_power_tot','LFP_power_tot');

% find common data
fw = sqrt(R.N(1)); %63;
hw = fw/2;
t_p = find(R.neuron_sample.t_ind{1}); % must has step size 1
if min(diff(t_p)) > 1
    warning('neuro sample should have a time step size = 1!s')
end
grid_is_common = R.grid.t_mid >= t_p(1) & R.grid.t_mid <= t_p(end);
t_common = R.grid.t_mid( grid_is_common );

% reformat spike data
mode = 'bayes';
switch mode
    case 'bayes'
        strut_tmp = R.grid.bayes;
    case 'quick'
        strut_tmp = R.grid.quick;
end
spike_x = strut_tmp.centre(1, grid_is_common );
spike_y = strut_tmp.centre(2, grid_is_common );
spike_w = strut_tmp.radius(:, grid_is_common );
spike_bf_log10 = strut_tmp.bayes_factor_ln(:, grid_is_common)/log(10); % change base
spike_h = strut_tmp.height(:, grid_is_common );
spike_tot = R.grid.num_spikes_win( grid_is_common );

% reformat ripple data
peak_common = peak( t_common - t_p(1) + 1, :  )'; %#ok<NODEF>
fit_common = ripple_fit( t_common - t_p(1) + 1  );
RPI_common = RPI(t_common - t_p(1) + 1);
ripple_power_tot_common = ripple_power_tot(t_common - t_p(1) + 1);
LFP_power_tot_common = LFP_power_tot(t_common - t_p(1) + 1);
peak_LFP_common = peak_LFP(t_common - t_p(1) + 1, :)'; %#ok<NODEF>
fit_goodness_common = ripple_fit_goodness( t_common - t_p(1) + 1  );
ripple_x = zeros(1,length(fit_common));
ripple_y = zeros(1,length(fit_common));
ripple_h = zeros(1,length(fit_common));
ripple_w = zeros(1,length(fit_common));
ripple_g = zeros(1,length(fit_common));
for i = 1:length(fit_common)
    ripple_x(1,i) = fit_common{i}.x_c - (round(fw/2) - peak_common(1,i)) - hw; %#ok<*SAGROW>
    if  ripple_x(1,i) < -hw
        ripple_x(1,i) =  ripple_x(1,i) + fw;
    elseif ripple_x(1,i) > hw
        ripple_x(1,i) = ripple_x(1,i) - fw;
    end
    
    ripple_y(1,i) = fit_common{i}.y_c - (round(fw/2) - peak_common(2,i)) - hw;
    if  ripple_y(1,i) < -hw
        ripple_y(1,i) =  ripple_y(1,i) + fw;
    elseif ripple_y(1,i) > hw
        ripple_y(1,i) = ripple_y(1,i) - fw;
    end
    ripple_h(1,i) = fit_common{i}.h;
    ripple_w(1,i) = fit_common{i}.sigma;
    ripple_g(1, i) = fit_goodness_common{i}.adjrsquare;
end

% discard bad fits
spike_bad = spike_bf_log10 < 2;
spike_x(spike_bad) = NaN;
spike_y(spike_bad) = NaN;
spike_h(spike_bad) = NaN;
spike_w(spike_bad) = NaN;

ripple_bad = ripple_g < 0.8; %(nanmean(ripple_g) - 2*nanstd(ripple_g));
ripple_x(ripple_bad) = NaN;
ripple_y(ripple_bad) = NaN;
ripple_h(ripple_bad) = NaN;
ripple_w(ripple_bad) = NaN;



% do position difference vs time lag analysis
mean_pos_diff_lag = [];
lag_steps =  -40 : 20;
lag_ms = lag_steps*dt*diff(t_common(1:2));
for i = 1:length(lag_steps) % positive lag means A lags behind B
    len = length(spike_x);
    lag_step_tmp = lag_steps(i);
    if lag_step_tmp < 0
        indA = 1:len+lag_step_tmp;
        indB = 1-lag_step_tmp:len;
    else
        indA = 1+lag_step_tmp:len;
        indB = 1:len-lag_step_tmp;
    end
    x_diff = min([abs(spike_x(indA)- ripple_x(indB)); 63-abs(spike_x(indA)- ripple_x(indB))]);
    y_diff = min([abs(spike_y(indA)- ripple_y(indB)); 63-abs(spike_y(indA)- ripple_y(indB))]);
    pos_diff = sqrt(x_diff.^2 + y_diff.^2);
    mean_pos_diff_lag = [mean_pos_diff_lag nanmean(pos_diff)]; %#ok<AGROW>
end

% find the optimal time lag
[~, min_pos_diff_ind] = min(mean_pos_diff_lag);
lag_ms_opt = lag_ms(min_pos_diff_ind);
lag_step_tmp = lag_steps(min_pos_diff_ind);
if lag_step_tmp < 0
    indA = 1:len+lag_step_tmp;
    indB = 1-lag_step_tmp:len;
else
    indA = 1+lag_step_tmp:len;
    indB = 1:len-lag_step_tmp;
end
x_diff = min([abs(spike_x(indA)- ripple_x(indB)); 63-abs(spike_x(indA)- ripple_x(indB))]);
y_diff = min([abs(spike_y(indA)- ripple_y(indB)); 63-abs(spike_y(indA)- ripple_y(indB))]);
pos_diff_opt = sqrt(x_diff.^2 + y_diff.^2);


% % collect spike hist according to raw ripple peak position
% n_bins = 5;
% t_tmp = t_common(indA);
% raw_h = peak_common(3,indB);
% raw_x = round(peak_common(1,indB));
% raw_y = round(peak_common(2,indB));
% spike_hist_t_range = -100:100; % steps
% spike_ripple_h_bin = linspace(min(raw_h), max(raw_h), n_bins+1);
% [x_grid_0, y_grid_0] = meshgrid(1:fw, 1:fw);
% spike_hist_acc_c = zeros(1,n_bins);
% spike_hist_acc = zeros(fw^2, length(spike_hist_t_range), n_bins);
%
% for i = 1:length(t_tmp)
%     if ~isnan( raw_x(i) ) && ~isnan( raw_h(i) )
%         t = t_tmp(i);
%         t_range_tmp = t+spike_hist_t_range;
%         if min(t_range_tmp) > 1 && max(t_range_tmp) <= R.step_tot
%             [ind_neu, ind_t] = find(R.spike_hist{1}(:, t_range_tmp));
%             if ~isempty(ind_neu)
%                 x_grid = circshift(x_grid_0, [round(fw/2-raw_x(i))  round(fw/2-raw_y(i))]);
%                 y_grid = circshift(y_grid_0, [round(fw/2-raw_x(i))  round(fw/2-raw_y(i))]);
%                 ind_neu_shifted = sub2ind([fw fw], x_grid(ind_neu), y_grid(ind_neu));
%                 bin_ind_tmp = find(raw_h(i) >= spike_ripple_h_bin(1:end-1) & raw_h(i) <= spike_ripple_h_bin(2:end), 1);
%                 spike_hist_tmp = full(sparse(ind_neu_shifted, ind_t, ones(size(ind_t)), fw^2, length(spike_hist_t_range)));
%                 spike_hist_acc(:,:,bin_ind_tmp) = spike_hist_acc(:,:,bin_ind_tmp) + spike_hist_tmp;
%                 spike_hist_acc_c(bin_ind_tmp) = spike_hist_acc_c(bin_ind_tmp) + 1;
%             end
%         end
%     end
% end

% collect spike hist according to raw ripple peak position
t_tmp = t_common(indA);
raw_h = peak_common(3,indB);
raw_x = round(peak_common(1,indB));
raw_y = round(peak_common(2,indB));
t_range_ms = 20; %ms;

spike_t_range_steps = round(t_range_ms/R.reduced.dt); %ms;

dt_conv = R.dt/R.reduced.dt;

[x_grid_0, y_grid_0] = meshgrid(1:fw, 1:fw);


[Lattice, ~] = lattice_nD(2, (fw-1)/2);
spike_sort_range = 6; % use 60
centre_neurons = sqrt( Lattice(:,1).^2 + Lattice(:,2).^2 ) <= spike_sort_range;

spike_count_sort = [];
spike_count_h = [];
for i = 1:length(t_tmp)
    if ~isnan( raw_x(i) ) && ~isnan( raw_h(i) )
        spike_count_sort_tmp = [];
        % spikes
        t = round(t_tmp(i)*dt_conv);
        t_range_tmp = (t-spike_t_range_steps):(t+spike_t_range_steps) ;
        if min(t_range_tmp) > 1 && max(t_range_tmp) <= R.reduced.step_tot
            
            [ind_neu, ind_t] = find(R.reduced.spike_hist{1}(:, t_range_tmp));
            
            x_grid = circshift(x_grid_0, [round(fw/2-raw_x(i))  round(fw/2-raw_y(i))]);
            y_grid = circshift(y_grid_0, [round(fw/2-raw_x(i))  round(fw/2-raw_y(i))]);
            ind_neu_shifted = sub2ind([fw fw], x_grid(ind_neu), y_grid(ind_neu));
            % plot(Lattice(ind_neu_shifted,1),Lattice(ind_neu_shifted,2),'o')
            spike_hist_tmp = full(sparse(ind_neu_shifted, ind_t, ones(size(ind_t)), fw^2, length(t_range_tmp)));
            spike_count_sort_tmp = sum(spike_hist_tmp(centre_neurons, :));
        end
        
        if ~isempty(spike_count_sort_tmp) %&& ~isempty(hil_phase_tmp)
            spike_count_sort = [spike_count_sort; spike_count_sort_tmp(:)']; %#ok<AGROW>
            spike_count_h = [spike_count_h; raw_h(i)]; %#ok<AGROW>
        end
    end
    
    
    
end


% output results
R.grid_SWR.spike_x = spike_x;
R.grid_SWR.spike_y = spike_y;
R.grid_SWR.spike_h = spike_h;
R.grid_SWR.spike_bf_log10 = spike_bf_log10;
R.grid_SWR.spike_w = spike_w;
R.grid_SWR.spike_tot = spike_tot;
R.grid_SWR.spike_ind = indA;

R.grid_SWR.ripple_x = ripple_x;
R.grid_SWR.ripple_y = ripple_y;
R.grid_SWR.ripple_h = ripple_h;
R.grid_SWR.ripple_g = ripple_g;
R.grid_SWR.ripple_w = ripple_w;
R.grid_SWR.ripple_ind = indB;

R.grid_SWR.ripple_power_index = RPI_common;
R.grid_SWR.ripple_power_tot = ripple_power_tot_common;
R.grid_SWR.LFP_power_tot = LFP_power_tot_common;
R.grid_SWR.LFP_h = peak_LFP_common(3,:);

R.grid_SWR.lag_ms = lag_ms;
R.grid_SWR.mean_pos_diff_lag = mean_pos_diff_lag;
R.grid_SWR.lag_ms_opt = lag_ms_opt;
R.grid_SWR.pos_diff_opt = pos_diff_opt;

% R.grid_SWR.spike_img_acc = spike_img_acc;
% R.grid_SWR.spike_img_acc_c = spike_img_acc_c;
% R.grid_SWR.spike_hist_acc = spike_hist_acc;
% R.grid_SWR.spike_hist_acc_c = spike_hist_acc_c;
% R.grid_SWR.spike_hist_t_range = spike_hist_t_range;
% R.grid_SWR.spike_ripple_h_bin = spike_ripple_h_bin;
R.grid_SWR.spike_count_sort = spike_count_sort;
R.grid_SWR.spike_count_h = spike_count_h;
end


