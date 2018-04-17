function [ R ] = get_grid_firing_centre( R, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% If you want to change mode,eg:[R] = get_grid_firing_centre(R,'mode','bayesian')

% parameters
win_len = 50; % window length in time steps
win_gap = 10; % window gap
win_min_rate_Hz = 0.5;

seg = 1:R.step_tot;
mode = 'quick'; % quick or bayesian
for i = 1:length(varargin)/2
    var_name = varargin{2*i-1};
    var_value = varargin{2*i};
     if isnumeric(var_value)
        eval([var_name, '=', num2str(var_value), ';']);
     else
         eval([var_name, '=''', var_value, ''';']);
     end
end
mode

spikes_win_min = win_min_rate_Hz*(R.dt*0.001)*win_len*R.N(1);
[ t_mid_full, ind_ab_full,  num_spikes_win_full, ~ ] = window_spike_hist_compressed( R, win_len, win_gap );


t_seg = t_mid_full > min(seg) &  t_mid_full <= max(seg);
t_mid = t_mid_full(t_seg);
ind_ab = ind_ab_full(:, t_seg);
num_spikes_win = num_spikes_win_full(t_seg);

    

% % spikes_win_min
% min_spike_requiremnt = num_spikes_win >= spikes_win_min;
% t_mid = t_mid( min_spike_requiremnt );
ind_a_vec = ind_ab(1,:);
ind_b_vec = ind_ab(2,:);
% % num_spikes_win = num_spikes_win( num_spikes_win >= spikes_win_min);


%%%% get window-ed mean and std for x/y position of firing neurons
hw = (R.N(1)^0.5 - 1)/2;
fw = 2*hw+1;


if mod(hw, 1) ~= 0
    % warning('Not a square grid')
else
    
    [Lattice, ~] = lattice_nD(2, hw);
    
    x_pos = Lattice(R.spike_hist_compressed{1}, 1);
    y_pos = Lattice(R.spike_hist_compressed{1}, 2);
    x_mean =  [];
    y_mean =  [];
    width =  [];
    mlh =  [];
    height = [];
    bayes_factor_ln = [];
    
    for j = 1:length(ind_a_vec)
        if mod(j*10,round(length(ind_a_vec)/10)*10) == 0
            fprintf('%d...', 10 - j*10 / (round(length(ind_a_vec)/10)*10));
        end
        % j/length(ind_a_vec)
        ind_range_tmp = ind_a_vec(j):ind_b_vec(j);
        x_pos_tmp = x_pos(ind_range_tmp);
        y_pos_tmp = y_pos(ind_range_tmp);
        
        if length(x_pos_tmp) < spikes_win_min
            x_mean_tmp = NaN;
            y_mean_tmp = NaN;
            width_tmp = NaN;
            mlh_tmp = NaN;
            height_tmp = NaN;
            bayes_factor_tmp = NaN;
        else
            [ x_mean_tmp, y_mean_tmp, width_tmp, mlh_tmp, height_tmp, bayes_factor_tmp ] = fit_bayesian_bump_2_spikes_circular(x_pos_tmp,y_pos_tmp, fw, mode);
        end
        
        x_mean =  [x_mean x_mean_tmp]; %#ok<AGROW>
        y_mean =  [y_mean y_mean_tmp];%#ok<AGROW>
        width =  [width width_tmp];%#ok<AGROW>
        mlh =  [mlh mlh_tmp];%#ok<AGROW>
        height = [height height_tmp]; %#ok<AGROW>
        bayes_factor_ln = [bayes_factor_ln bayes_factor_tmp]; %#ok<AGROW>
    end
    
end

% deal with periodic boundary
x_s = [-fw 0 fw];
y_s = [-fw 0 fw];
[x_s_grid, y_s_grid] = meshgrid(x_s, y_s);
x_shift = x_s_grid(:);
y_shift = y_s_grid(:);
no_shifts = length(y_shift);

% raw jump_dir
x_diff_full = repmat(x_mean(2:end), no_shifts, 1);
y_diff_full = repmat(y_mean(2:end), no_shifts, 1);
for s = 1:no_shifts
    x_diff_full(s,:) =  x_diff_full(s,:) - x_shift(s) - x_mean(1:end-1);
    y_diff_full(s,:) =  y_diff_full(s,:) - y_shift(s) - y_mean(1:end-1) ;
end
[jump_dist_raw, min_J] = min(sqrt(x_diff_full.^2 + y_diff_full.^2));
jump_dir_raw = zeros(size(jump_dist_raw));
for i = 1:length(min_J)
    jump_dir_raw(i) = atan2(y_diff_full(min_J(i),i),  x_diff_full(min_J(i),i));
end

% calculate mean jerk
dt_mid = t_mid_full(2) - t_mid_full(1);
x_tmp = cos(jump_dir_raw).*jump_dist_raw;
y_tmp = sin(jump_dir_raw).*jump_dist_raw;
x1 = diff(x_tmp)/dt_mid;
x2 = diff(x1)/dt_mid;
x3 = diff(x2)/dt_mid;
y1 = diff(y_tmp)/dt_mid;
y2 = diff(y1)/dt_mid;
y3 = diff(y2)/dt_mid;
jerk_mean = nanmean(sqrt(x3.^2 + y3.^2));



% output results
R.grid.win_len = win_len;
R.grid.win_gap = win_gap; % window gap
R.grid.spikes_win_min =  spikes_win_min;
R.grid.win_min_rate_Hz = win_min_rate_Hz;
R.grid.num_spikes_win = num_spikes_win_full;
R.grid.t_mid = t_mid_full;
R.grid.ind_ab = ind_ab;
switch lower(mode)
    case 'quick'
        R.grid.quick.radius = width;
        R.grid.quick.centre = [x_mean; y_mean];
        R.grid.quick.jump_dist = jump_dist_raw;
        R.grid.quick.jump_dir = jump_dir_raw;
        R.grid.quick.mlh = mlh;
        R.grid.quick.height = height;
        R.grid.quick.jerk_mean = jerk_mean;
    case 'bayesian'
        R.grid.bayes.radius = width;
        R.grid.bayes.centre = [x_mean; y_mean];
        R.grid.bayes.jump_dist = jump_dist_raw;
        R.grid.bayes.jump_dir = jump_dir_raw;
        R.grid.bayes.mlh = mlh;
        R.grid.bayes.height = height;
        R.grid.bayes.jerk_mean = jerk_mean;
        R.grid.bayes.bayes_factor_ln = bayes_factor_ln;
end

if ~isfield(R, 'grid_sub') 
    R_sub.grid_sub = []; % to stop recursion
    % do the substitute study
    n_s = sum(R.spike_hist{1},2);
    i_s = [];
    t_s = [];
    for i = 1:R.N(1)
        i_s = [i_s i*ones(1,n_s(i))]; %#ok<*AGROW>
        t_s = [t_s randperm(R.step_tot, n_s(i))];
    end
    sh = sparse( i_s, t_s, ones(size(i_s)),  R.N(1), R.step_tot);
    R_sub.num_spikes{1} = sum(sh, 1);
    R_sub.N = R.N;
    R_sub.dt = R.dt;
    R_sub.step_tot = R.step_tot;
    [sc,~] = find(sh);
    R_sub.spike_hist_compressed{1} = sc;
    % 
    [ R_sub ] = get_grid_firing_centre( R_sub, varargin );
    R.grid_sub = R_sub.grid;
end

end

% function 
