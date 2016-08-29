function [ R ] = get_grid_firing_centre( R, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% parameters
win_len = 100; % window length in time steps
win_gap = 50; % window gap
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


t_seg = t_mid_full > min(seg) &  t_mid_full < max(seg);
t_mid = t_mid_full(t_seg);
ind_ab = ind_ab_full(:, t_seg);
num_spikes_win = num_spikes_win_full(t_seg);

    

% spikes_win_min
t_mid = t_mid( num_spikes_win >= spikes_win_min);
ind_a_vec = ind_ab(1,num_spikes_win >= spikes_win_min);
ind_b_vec = ind_ab(2,num_spikes_win >= spikes_win_min);
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
        
        
        % add bayesian stuff here
        [ x_mean_tmp, y_mean_tmp, width_tmp, mlh_tmp ] = fit_bayesian_bump_2_spikes(x_pos_tmp,y_pos_tmp, fw, mode);
        x_mean =  [x_mean x_mean_tmp]; %#ok<AGROW>
        y_mean =  [y_mean y_mean_tmp];%#ok<AGROW>
        width =  [width width_tmp];%#ok<AGROW>
        mlh =  [mlh mlh_tmp];%#ok<AGROW>
        
    end
    
end

%thre
c_legit = (mlh >= mean(mlh) - 2*std(mlh)); % this is a bit arbitarys!!!

t_mid_chosen = t_mid(c_legit);
width_chosen = width(c_legit);
x_mean_chosen = x_mean(c_legit);
y_mean_chosen = y_mean(c_legit);

% jump_size (take care of the periodic boundary condition)
x_diff_full = repmat(x_mean_chosen(2:end), 9, 1);
y_diff_full = repmat(y_mean_chosen(2:end), 9, 1);
x_mean_shift2 = [-fw 0  fw -fw 0 fw -fw  0   fw];
y_mean_shift2 = [ fw fw fw   0 0  0 -fw -fw -fw];
for s = 1:9
    x_diff_full(s,:) =  x_diff_full(s,:) - x_mean_shift2(s) - x_mean_chosen(1:end-1);
    y_diff_full(s,:) =  y_diff_full(s,:) - y_mean_shift2(s) - y_mean_chosen(1:end-1) ;
end
jump_dist = min(sqrt(x_diff_full.^2 + y_diff_full.^2));
% jump_duration
jump_duration = diff(t_mid_chosen);

% output results

R.grid.raw.win_len = win_len;
R.grid.raw.win_gap = win_gap; % window gap
R.grid.raw.win_min_rate_Hz = win_min_rate_Hz;
R.grid.raw.num_spikes_win = num_spikes_win_full;
R.grid.raw.t_mid = t_mid_full;
R.grid.raw.ind_ab = ind_ab_full;

% R.grid.raw.t_ab = t_ab_vec_full; 

R.grid.mode = mode;
R.grid.t_mid = t_mid_chosen;
R.grid.radius = width_chosen;
R.grid.centre = [x_mean_chosen; y_mean_chosen];
R.grid.mlh = mlh;
R.grid.jump_dist = jump_dist;
R.grid.jump_du = jump_duration;

end


