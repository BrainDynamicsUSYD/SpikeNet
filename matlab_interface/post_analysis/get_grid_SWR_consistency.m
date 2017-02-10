function R = get_grid_SWR_consistency(R)

% Load data
dt = R.dt;
stamp = R.stamp;
samp_file = [stamp(1:end-3) '0_neurosamp'];
load(['./', samp_file], 'peak','ripple_fit', 'ripple_fit_goodness');

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
        
% reformat ripple data
peak_common = peak( t_common - t_p(1) + 1, :  )';
fit_common = ripple_fit( t_common - t_p(1) + 1  );
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

ripple_bad = ripple_g < (nanmean(ripple_g) - 2*nanstd(ripple_g));
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


% collect spike image according to ripple_h_bin
n_bins = 5;
ripple_h_tmp = ripple_h(indB);
[Lattice, ~] = lattice_nD(2, (fw-1)/2);
t_tmp = t_common(indA);
spike_x_tmp = round(spike_x(indA)+fw/2);
spike_y_tmp = round(spike_y(indA)+fw/2);
spike_img_t_range = -25:25; % steps
spike_img_ripple_h_bin = linspace(min(ripple_h_tmp), max(ripple_h_tmp), n_bins+1);
spike_img_acc_c = zeros(1,n_bins);
spike_img_acc = zeros(fw, fw, n_bins);
for i = 1:length(t_tmp)
    if ~isnan( spike_x_tmp(i) ) && ~isnan(ripple_h_tmp(i))
        t = t_tmp(i);
        t_range_tmp = t+spike_img_t_range;
        if min(t_range_tmp) > 1 && max(t_range_tmp) <= R.step_tot
            [ind_neu, ind_t] = find(R.spike_hist{1}(:, t_range_tmp));
            ind_neu_unique = unique(ind_neu);
            i_tmp = Lattice(ind_neu_unique,1)+32;
            j_tmp = Lattice(ind_neu_unique,2)+32;
            ind_neu_occ = histc(ind_neu, ind_neu_unique);
            img_tmp = full(sparse(i_tmp, j_tmp, ind_neu_occ, fw,fw));
            % shift to the center
            img_tmp = circshift(img_tmp, [round(fw/2)-spike_x_tmp(i)  round(fw/2)-spike_y_tmp(i)]);
            bin_ind_tmp = find(ripple_h_tmp(i) >= spike_img_ripple_h_bin(1:end-1) & ripple_h_tmp(i) <= spike_img_ripple_h_bin(2:end));
            spike_img_acc(:,:,bin_ind_tmp) = spike_img_acc(:,:,bin_ind_tmp) + img_tmp;
            spike_img_acc_c(bin_ind_tmp) = spike_img_acc_c(bin_ind_tmp) + 1;
        end
    end
end


% output results
R.grid_SWR.spike_x = spike_x;
R.grid_SWR.spike_y = spike_y;
R.grid_SWR.spike_h = spike_h;
R.grid_SWR.spike_bf_log10 = spike_bf_log10;
R.grid_SWR.spike_w = spike_w;
R.grid_SWR.spike_ind = indA;

R.grid_SWR.ripple_x = ripple_x;
R.grid_SWR.ripple_y = ripple_y;
R.grid_SWR.ripple_h = ripple_h;
R.grid_SWR.ripple_g = ripple_g;
R.grid_SWR.ripple_w = ripple_w;
R.grid_SWR.ripple_ind = indB;

R.grid_SWR.lag_ms = lag_ms;
R.grid_SWR.mean_pos_diff_lag = mean_pos_diff_lag;
R.grid_SWR.lag_ms_opt = lag_ms_opt;
R.grid_SWR.pos_diff_opt = pos_diff_opt;

R.grid_SWR.spike_img_acc = spike_img_acc;
R.grid_SWR.spike_img_acc_c = spike_img_acc_c;
R.grid_SWR.spike_img_t_range = spike_img_t_range; % steps
R.grid_SWR.spike_img_ripple_h_bin = spike_img_ripple_h_bin;
end


