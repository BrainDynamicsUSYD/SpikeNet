function [ t_mid, ind_ab,  num_spikes_win, t_ab_vec ] = window_spike_hist_compressed( R, win_len, win_gap, pop_ind )
% [ t_mid, ind_ab,  num_spikes_win ] = window_spike_hist_compressed( R, win_len, win_gap, pop_ind )
%This function returns ind_ab = [ind_a; ind_b] so that
%spike_hist_compressed{pop_ind}(ind_a(i):ind_b(i)) gives the indices of the
%neurons that spiked during the i-th window as defined by win_len and
%win_gap (unit: steps). t_mid is the centre of the windows.


if nargin == 3
    pop_ind = 1;
end


%%%% get window-ed index range for spike_hist_compressed
t_a_vec = [];
t_b_vec = [];
ind_a_vec = [];
ind_b_vec = [];
num_spikes_win = [];

t_a = 1;
t_b = win_len;
ind_a = 1;
while t_b  < R.step_tot
    t_a_vec = [t_a_vec t_a];%#ok<AGROW>
    t_b_vec = [t_b_vec t_b];%#ok<AGROW>
    
    num_spikes_win_tmp = sum(R.num_spikes{pop_ind}(t_a:t_b));
    num_spikes_win = [num_spikes_win num_spikes_win_tmp];%#ok<AGROW>
    ind_b = ind_a + num_spikes_win_tmp - 1;
    
    ind_a_vec = [ind_a_vec ind_a];%#ok<AGROW>
    ind_b_vec = [ind_b_vec ind_b];%#ok<AGROW>
    
    num_spikes_win_gap = sum(R.num_spikes{pop_ind}(t_a:(t_a+win_gap-1)));
    ind_a = ind_a + num_spikes_win_gap;
    
    t_a = t_a + win_gap;
    t_b = t_b + win_gap;
end
t_ab_vec = [t_a_vec(:) t_b_vec(:)];
t_mid = round( (t_a_vec+t_b_vec)/2 );
ind_ab = [ind_a_vec; ind_b_vec];

% % The above code has been tested using the following code
% t_s = 1234;
% nnz(R.spike_hist{pop_ind}(:,t_a_vec(t_s):t_b_vec(t_s)))
% num_spikes_win(t_s)
% nnz(R.spike_hist_compressed{pop_ind}(ind_a_vec(t_s):ind_b_vec(t_s)))


end

