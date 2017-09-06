function R = get_triplet_sequence_no_latency(R, varargin)
fprintf('Getting triplet results...');

sh0 = R.reduced.spike_hist{1};
dt = R.reduced.dt; % 1 ms

% up-state duration
up_du = 100;

% 2D filter
bin = 3.2; % sec
gauss_sigma = 10; % sec
gauss_sigma_bin = gauss_sigma/bin;
max_t_diff = 300;
ed = -max_t_diff:bin:max_t_diff;
C = {ed, ed};

% sampling area of silicon microelectrode
n_trial = 50; %20;
n_trip = 50;
% n_circle = 1000;
n_sample = 50;
hw_sample = 5;

for i = 1:length(varargin)/2
    var_name = varargin{2*i-1};
    var_value = varargin{2*i};
     if isnumeric(var_value)
        eval([var_name, '=', num2str(var_value), ';']);
     else
         eval([var_name, '=''', var_value, ''';']);
     end
end


L = lattice_nD(2, 31);
% I_sample = find(  L(:,1).^2+L(:,2).^2 <= n_circle/pi ); % circular method
I_sample = find(  L(:,1).^2 <= hw_sample^2 &  L(:,2).^2 <= 25^2 ); % rectangular method
sh0 = sh0(I_sample,:);



% sampling neurons
sh_c = cell(1,n_trial);
sh_ind = cell(1,n_trial); 
for jj = 1:n_trial % need this to get a wide range of latency
    ind_rand = randperm(length(I_sample),n_sample);
    sh_c{jj} = sh0(ind_rand, :);
    sh_ind{jj} = I_sample(ind_rand);
end


up_onset_c = cell(1,n_trial);
for jj = 1:n_trial % need this to get a wide range of latency
    
    sh = sh_c{jj};
    % up and down state detection
    [up_onset, up_offset] = get_up_and_down( sh, dt,'up_du', up_du );
    n_on = length(up_onset);

    up_onset_c{jj} = up_onset;
   
    
end


ind_trip_c = cell(1, n_trial);
N_trip_c = cell(1, n_trial);
trip_peak_c = cell(1, n_trial);
N_trip_shuffle_c = cell(1, n_trial);
trip_peak_shuffle_c = cell(1, n_trial);
t_from_onset_shuffle_c =  cell(1, n_trial);
t_from_onset_c{jj} =  cell(1, n_trial);
for jj = 1:n_trial
    jj/n_trial
    sh = sh_c{jj};
    up_onset = up_onset_c{jj};
    trip_peak = zeros(n_trip, 3);
    ind_trip_mat = zeros(n_trip, 3);
    t_from_onset = cell(1,n_trip);
    for i = 1:n_trip
%         ah(i) = subaxis(6,6,i,'PR',0.01);
        ind_trip = randperm(n_sample, 3);
        
        ind_trip_mat(i, :) = ind_trip(:)';
        sp = sh(ind_trip, :);
        
        %         [t1, t2, t3] = meshgrid( find(sp(1,:)), find(sp(2,:)), find(sp(3,:)) );
        [t1, t2, t3] = get_t123(sp, max_t_diff);
        ta = t1-t3;
        tb = t2-t3;
        [N_trip,~] = hist3([ta(:), tb(:)],'Edges',{ed, ed});

        N_trip = N_trip / (sum(sp(3,:)) * (bin/1000)^2 );
        N_trip = imgaussfilt(N_trip,gauss_sigma_bin) ;

        %     imagesc(C{1},C{2}, N)
        %     axis([-150 150 -150 150])
        
        [n_max,n_ind] = max(N_trip(:));
        [n_i,n_j] = ind2sub(size(N_trip),n_ind);
        trip_peak(i,:) = [C{1}(n_i), C{2}(n_j), n_max];
        
        % time from onset 
        
        is_trip = (abs(ta - C{1}(n_i)) < 10 & abs(tb - C{2}(n_j)) < 10);
        t1_is_trip = t1(is_trip);
        t2_is_trip = t2(is_trip);
        t3_is_trip = t3(is_trip);
        t_trip = min( [t1_is_trip(:)'; t2_is_trip(:)'; t3_is_trip(:)']  );
        t_from_onset_tmp =  t_trip(:)' -  up_onset(:);
        t_from_onset_tmp(t_from_onset_tmp < 0) = NaN;
        t_from_onset_tmp = min(t_from_onset_tmp);
        t_from_onset{i} = t_from_onset_tmp(:)'; 
    end
    % set(gca,'CLIM',[0 320])
    % colormap('jet')
    
    % shuffled triplet
    trip_peak_shuffle = zeros(n_trip, 3);
    sh_shuffle = raster_marginals_shuffling(sh);
    t_from_onset_shuffle = cell(1,n_trip);
    for i = 1:n_trip
%         ah(i) = subaxis(6,6,i,'PR',0.01);
        ind_trip = ind_trip_mat(i,:);
        sp = sh_shuffle(ind_trip, :);
        
        %         [t1, t2, t3] = meshgrid( find(sp(1,:)), find(sp(2,:)), find(sp(3,:)) );
        [t1, t2, t3] = get_t123(sp, max_t_diff);
        ta = t1-t3;
        tb = t2-t3;
        [N_trip_s,~] = hist3([ta(:), tb(:)],'Edges',{ed, ed});

        N_trip_s = N_trip_s / (sum(sp(3,:)) * (bin/1000)^2 );
        N_trip_s = imgaussfilt(N_trip_s,gauss_sigma_bin) ;
        
        
        %     imagesc(C{1},C{2}, N)
        %     axis([-150 150 -150 150])
        
        [n_max,n_ind] = max(N_trip_s(:));
        [n_i,n_j] = ind2sub(size(N_trip_s),n_ind);
        trip_peak_shuffle(i,:) = [C{1}(n_i), C{2}(n_j),n_max];
        
        
         % time from onset 
        is_trip = (abs(ta - C{1}(n_i)) < 10 & abs(tb - C{2}(n_j)) < 10);
        t1_is_trip = t1(is_trip);
        t2_is_trip = t2(is_trip);
        t3_is_trip = t3(is_trip);
        t_trip = min( [t1_is_trip(:)'; t2_is_trip(:)'; t3_is_trip(:)']  );
        t_from_onset_tmp =  t_trip(:)' -  up_onset(:);
        t_from_onset_tmp(t_from_onset_tmp < 0) = NaN;
        t_from_onset_tmp = min(t_from_onset_tmp);
        t_from_onset_shuffle{i} = t_from_onset_tmp(:)'; 
    end
    % set(gca,'CLIM',[0 320])
    % colormap('jet')

    ind_trip_c{jj} = ind_trip_mat;
    N_trip_c{jj} = N_trip; % only record the last one to save some memory
    trip_peak_c{jj} = trip_peak;
    N_trip_shuffle_c{jj} = N_trip_s; %  only record the last one to save some memory
    trip_peak_shuffle_c{jj} = trip_peak_shuffle;
    t_from_onset_shuffle_c{jj} =  t_from_onset_shuffle;
    t_from_onset_c{jj} =  t_from_onset;
end


R.triplet.sh_ind_c = sh_ind;
R.triplet.ind_trip_c = ind_trip_c;
R.triplet.N_trip_c = N_trip_c;
R.triplet.N_trip_bin = C;
R.triplet.trip_peak_c = trip_peak_c;
R.triplet.N_trip_shuffle_c = N_trip_shuffle_c;
R.triplet.trip_peak_shuffle_c = trip_peak_shuffle_c;
R.triplet.t_from_onset_shuffle_c =  t_from_onset_shuffle_c;
R.triplet.t_from_onset_c =  t_from_onset_c;
disp('Done.');

end


function [t1, t2, t3] = get_t123(sp, max_t_diff)
t1 = [];
t2 = [];
t3 = [];

t3_all = find(sp(3,:));
for i = 1:length(t3_all)
    t3_tmp = t3_all(i);
    if t3_tmp - max_t_diff > 0 && t3_tmp + max_t_diff <= length(sp(3,:))
        t_range = (t3_tmp - max_t_diff):(t3_tmp + max_t_diff);
        
        
        [t1_tmp, t2_tmp] = meshgrid( find(sp(1,t_range)), find(sp(2,t_range)) );

        t1 = [t1 t1_tmp(:)']; %#ok<*AGROW>
        t2 = [t2 t2_tmp(:)'];
        t3 = [t3 t3_tmp*ones(size(t2_tmp(:)'))];
    end
    
end

end
