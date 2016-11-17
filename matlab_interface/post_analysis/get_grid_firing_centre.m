function [ R ] = get_grid_firing_centre( R, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


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

    

% spikes_win_min
min_spike_requiremnt = num_spikes_win >= spikes_win_min;
t_mid = t_mid( min_spike_requiremnt );
ind_a_vec = ind_ab(1,min_spike_requiremnt);
ind_b_vec = ind_ab(2,min_spike_requiremnt);
% num_spikes_win = num_spikes_win( num_spikes_win >= spikes_win_min);


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
    
    for j = 1:length(ind_a_vec)
        if mod(j*10,round(length(ind_a_vec)/10)*10) == 0
            fprintf('%d...', 10 - j*10 / (round(length(ind_a_vec)/10)*10));
        end
        % j/length(ind_a_vec)
        ind_range_tmp = ind_a_vec(j):ind_b_vec(j);
        x_pos_tmp = x_pos(ind_range_tmp);
        y_pos_tmp = y_pos(ind_range_tmp);
        
%         if j> 298
%             j
%         end
        % add bayesian stuff here
        [ x_mean_tmp, y_mean_tmp, width_tmp, mlh_tmp ] = fit_bayesian_bump_2_spikes(x_pos_tmp,y_pos_tmp, fw, mode);
        x_mean =  [x_mean x_mean_tmp]; %#ok<AGROW>
        y_mean =  [y_mean y_mean_tmp];%#ok<AGROW>
        width =  [width width_tmp];%#ok<AGROW>
        mlh =  [mlh mlh_tmp];%#ok<AGROW>
        
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
for s = 1:no_shifts;
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
jerk_mean = mean(sqrt(x3.^2 + y3.^2));

%threholding on MLH
c_legit = (mlh >= mean(mlh) - 2*std(mlh)); % this is a bit arbitarys!!!

t_mid_chosen = t_mid(c_legit);
width_chosen = width(c_legit);
x_mean_chosen = x_mean(c_legit);
y_mean_chosen = y_mean(c_legit);

% jump_size (take care of the periodic boundary condition)
x_diff_full = repmat(x_mean_chosen(2:end), no_shifts, 1);
y_diff_full = repmat(y_mean_chosen(2:end), no_shifts, 1);

for s = 1:no_shifts;
    x_diff_full(s,:) =  x_diff_full(s,:) - x_shift(s) - x_mean_chosen(1:end-1);
    y_diff_full(s,:) =  y_diff_full(s,:) - y_shift(s) - y_mean_chosen(1:end-1) ;
end
[jump_dist, min_J] = min(sqrt(x_diff_full.^2 + y_diff_full.^2));
jump_dir = zeros(size(jump_dist));
for i = 1:length(min_J)
    jump_dir(i) = atan2(y_diff_full(min_J(i),i),  x_diff_full(min_J(i),i));
end

% jump_duration
jump_duration = diff(t_mid_chosen);

% output results

R.grid.raw.win_len = win_len;
R.grid.raw.win_gap = win_gap; % window gap
R.grid.raw.min_spike_requiremnt = min_spike_requiremnt;
R.grid.raw.win_min_rate_Hz = win_min_rate_Hz;
R.grid.raw.num_spikes_win = num_spikes_win_full;
R.grid.raw.t_mid = t_mid_full;
R.grid.raw.ind_ab = ind_ab_full;
R.grid.raw.radius = width;
R.grid.raw.centre = [x_mean; y_mean];
R.grid.raw.jump_dist = jump_dist_raw;
R.grid.raw.jump_dir = jump_dir_raw;
R.grid.raw.mlh = mlh;
R.grid.raw.jerk_mean = jerk_mean;

R.grid.mode = mode;
R.grid.t_mid = t_mid_chosen;
R.grid.radius = width_chosen;
R.grid.centre = [x_mean_chosen; y_mean_chosen];
R.grid.jump_dist = jump_dist;
R.grid.jump_du = jump_duration;
R.grid.jump_dir = jump_dir;

end


