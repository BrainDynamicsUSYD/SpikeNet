function [R] = get_fano_factor(R, spatial_temporal_option, hw)
fprintf('\t Getting Fano factor...\n');

if nargin == 1
    spatial_temporal_option = 1;
end

dt = R.reduced.dt;
spike_hist = R.reduced.spike_hist{1}'; % make time columnwise
N = R.N;
step_tot = R.reduced.step_tot;


trans_ms = 500; % ms, trasient period to be discarded
win_gap_ratio = 2; % should be larger than 1
ker_gap_ratio = 2; % should be larger than 1


if spatial_temporal_option == 1
    n_sample = 1000;
    ind_sample = randperm(N(1), n_sample);
    win_len = [10:10:100 200:100:1000]; % ms
    % temporal fano factor
    fano_tmpo = [];
    for win_tmp = win_len
        win_step = round(win_tmp/dt); % in steps
        trans_step = round(trans_ms/dt);
        t_b = trans_step+win_step+1 : round(win_gap_ratio*win_step): step_tot;
        t_a = t_b - win_step + 1;
        SC_tmp = zeros(length(t_a),n_sample);
        for i = 1:length(t_a)
            SC_tmp(i,:) = full(sum(spike_hist(t_a(i):t_b(i),ind_sample)));
        end
        fano_tmpo = [fano_tmpo nanmean(var(SC_tmp)./mean(SC_tmp))]; %#ok<AGROW>
    end
    
    R.Analysis.fano_tmpo = fano_tmpo;
    R.Analysis.fano_tmpo_win = win_len;
    R.Analysis.fano_n_sample = n_sample;
    
elseif spatial_temporal_option == 2
    ker_len = 4:2:16;
    spike_hist_full = R.spike_hist{1}';
    % spatial fano factor
    if nargin < 3
        hw = R.ExplVar.hw;
    end
    w = hw*2 + 1;
    
    for k = 1:length(ker_len)
        ker_len_tmp = ker_len(k);
        
        ker_neu_ind = [];
        
        xy_b = 1+ker_len_tmp:round(ker_len_tmp*ker_gap_ratio):w;
        xy_a = xy_b - ker_len_tmp + 1;
        
        [xy_a_x, xy_a_y] = meshgrid(xy_a,xy_a);
        [xy_b_x, xy_b_y] = meshgrid(xy_b,xy_b);
        
        for i = 1:length(xy_a)
            for j = 1:length(xy_a)
                g = zeros(w,w);
                g(xy_a_x(i,j):xy_b_x(i,j), xy_a_y(i,j):xy_b_y(i,j)) = 1;
                ker_neu_ind = [ker_neu_ind; find(g(:))']; %#ok<AGROW>
                %                 % debug
                %                 g_tmp = zeros(w,w);
                %                 g_tmp(ker_neu_ind(end,:)) = 1;
                %                 imagesc(g_tmp);
                %                 pause(0.1)
            end
        end
        
        
        
        SC_tmp = [];
        for i = 1:length(ker_neu_ind(:,1))
            SC_tmp  = [SC_tmp; sum(spike_hist_full(:,ker_neu_ind(i,:)),2)']; %#ok<AGROW>
            % disp(((i-1)*length(t_a) + t)/(length(ker_neu_ind(:,1))*length(t_a)))
        end
        fano_spat(k) =  nanmean(var(SC_tmp)./mean(SC_tmp)); %#ok<AGROW>
    end
    
    R.Analysis.fano_spat = fano_spat;
    R.Analysis.fano_spat_ker = ker_len;
    
elseif spatial_temporal_option == 3
    win_len = 50:50:500; % ms
    ker_len = 4:2:16;
    % spatial-tmporal fano factor
    if nargin < 3
        hw = R.ExplVar.hw;
    end
    fw = hw*2 + 1;
    
    for k = 1:length(ker_len)
        ker_len_tmp = ker_len(k);
        
        ker_neu_ind = [];
        
        xy_b = 1+ker_len_tmp:round(ker_len_tmp*ker_gap_ratio):fw;
        xy_a = xy_b - ker_len_tmp;
        
        [xy_a_x, xy_a_y] = meshgrid(xy_a,xy_a);
        [xy_b_x, xy_b_y] = meshgrid(xy_b,xy_b);
        
        for i = 1:length(xy_a)
            for j = 1:length(xy_a)
                g = zeros(fw,fw);
                g(xy_a_x(i,j):xy_b_x(i,j), xy_a_y(i,j):xy_b_y(i,j)) = 1;
                ker_neu_ind = [ker_neu_ind; find(g(:))']; %#ok<AGROW>
                %                 % debug
                %                 g_tmp = zeros(w,w);
                %                 g_tmp(ker_neu_ind(end,:)) = 1;
                %                 imagesc(g_tmp);
                %                 pause(0.1)
            end
        end
        
        
        for w = 1:length(win_len)
            win_tmp = win_len(w);
            win_step = round(win_tmp/dt); % in steps
            trans_step = round(trans_ms/dt);
            t_b = trans_step+win_step+1 : round(win_gap_ratio*win_step): step_tot;
            t_a = t_b - win_step + 1;
            
            SC_tmp = [];
            length(ker_neu_ind(:,1))
            for i = 1:length(ker_neu_ind(:,1))
                for t = 1:length(t_a)
                    SC_tmp  = [SC_tmp; sum(spike_hist(t_a(t):t_b(t),ker_neu_ind(i,:)))]; %#ok<AGROW>
                    % disp(((i-1)*length(t_a) + t)/(length(ker_neu_ind(:,1))*length(t_a)))
                end
            end
            fano_spat_tmpo(k,w) =  mean(var(SC_tmp)./mean(SC_tmp)); %#ok<AGROW>
            
        end
        
    end
    
    R.Analysis.fano_spat_tmpo = fano_spat_tmpo;
    R.Analysis.fano_win = win_len;
    R.Analysis.fano_ker = ker_len;
    
end

end





