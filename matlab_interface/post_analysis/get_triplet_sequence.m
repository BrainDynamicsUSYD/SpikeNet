function R = get_triplet_sequence(R, varargin)
fprintf('Getting triplet results...');

sh0 = R.reduced.spike_hist{1};
dt = R.reduced.dt; % 1 ms

% up-state duration
up_du = 100;

% 2D filter
bin = 3.2; % sec
gauss_sigma = 10; % sec
gauss_sigma_bin = gauss_sigma/bin;
ed = -300:bin:300;

% 1D filter
sigma = 10;
sz = 60;    % length of gaussFilter vector
x = linspace(-sz / 2, sz / 2, sz);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize


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

PETH_1st_c = cell(1,n_trial);
PETH_2nd_c = cell(1,n_trial);
latency_1st_c =  cell(1,n_trial);
latency_2nd_c =  cell(1,n_trial);
PETH_c = cell(1,n_trial);
latency_c = cell(1,n_trial);
PETH_unique_c = cell(1,n_trial);
up_onset_c = cell(1,n_trial);
up_offset_c = cell(1,n_trial);
for jj = 1:n_trial % need this to get a wide range of latency
    
    sh = sh_c{jj};
    % up and down state detection
    [up_onset, up_offset] = get_up_and_down( sh, dt,'up_du', up_du );
    n_on = length(up_onset);
    
    % PETH and latency
    PETH_1st = zeros(n_sample, up_du);
    PETH_2nd = zeros(n_sample, up_du);
    latency_1st = zeros(1,n_sample);
    latency_2nd = zeros(1,n_sample);
    for i = 1:n_sample
        for t = 1:round(n_on/2)
            a = sh(i,up_onset(t):up_onset(t)+up_du-1);
            PETH_1st(i,:) = PETH_1st(i,:) + a;
        end
        
        PETH_1st(i,:) = conv (PETH_1st(i,:), gaussFilter, 'same');
        latency_1st(i) = sum(PETH_1st(i,:).*(1:up_du))/sum(PETH_1st(i,:));  %#ok<*SAGROW>
        for t =  round(n_on/2)+1:n_on
            a = sh(i,up_onset(t):up_onset(t)+up_du-1);
            PETH_2nd(i,:) = PETH_2nd(i,:) + a;
        end
        
        PETH_2nd(i,:) = conv (PETH_2nd(i,:), gaussFilter, 'same');
        latency_2nd(i) = sum(PETH_2nd(i,:).*(1:up_du))/sum(PETH_2nd(i,:));
    end
    
    
    % PETH and latency
    PETH = zeros(n_sample, up_du);
    latency = zeros(1,n_sample);
    for i = 1:n_sample
        for t = 1:n_on
            a = sh(i,up_onset(t):up_onset(t)+up_du-1);
            PETH(i,:) = PETH(i,:) + a;
        end
        
        PETH(i,:) = conv (PETH(i,:), gaussFilter, 'same');
        latency(i) = sum(PETH(i,:).*(1:up_du))/sum(PETH(i,:));  %#ok<*SAGROW>
    end
    
    % PETH uniqueness
    PETH_unique = zeros(1,n_sample);
    for i = 1:n_sample
        D_i = zeros(1,n_sample);
        ai = PETH_1st(i,:)/sum(PETH_1st(i,:));
        for j = 1:n_sample
            aj = PETH_2nd(j,:)/sum(PETH_2nd(j,:));
            D_i(j) = sqrt(sum((ai(:) - aj(:)).^2));
        end
        PETH_unique(i) = sum(D_i >= D_i(i)) / n_sample;
    end
    
    up_onset_c{jj} = up_onset;
    up_offset_c{jj} = up_offset;
    PETH_1st_c{jj} = PETH_1st;
    PETH_2nd_c{jj} = PETH_2nd;
    latency_1st_c{jj} = latency_1st;
    latency_2nd_c{jj} = latency_2nd;
    PETH_c{jj} = PETH;
    latency_c{jj} = latency;
    PETH_unique_c{jj} =  PETH_unique;
    
end

% figure(1)
% hold on;
% plot(cell2mat(latency_1st_c(:)'), cell2mat(latency_2nd_c(:)'),'.')
% [C,P] = corrcoef(cell2mat(latency_1st_c(:)'),cell2mat(latency_2nd_c(:)'),'rows','complete') %#ok<*NOPTS>
% 
% mean(cell2mat(PETH_unique_c(:)'))
% mode(cell2mat(PETH_unique_c(:)'))

% triplet
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
        
        [t1, t2, t3] = meshgrid( find(sp(1,:)), find(sp(2,:)), find(sp(3,:)) );
        ta = t1-t3;
        tb = t2-t3;
        
        [N_trip,C] = hist3([ta(:), tb(:)],'Edges',{ed, ed});
        N_trip = N_trip / (sum(sp(3,:)) * (bin/1000)^2 );
        N_trip = imgaussfilt(N_trip,gauss_sigma_bin);
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
        t_from_onset{i} = t_from_onset_tmp(:)'; %#ok<AGROW>
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
        
        [t1, t2, t3] = meshgrid( find(sp(1,:)), find(sp(2,:)), find(sp(3,:)) );
        ta = t1-t3;
        tb = t2-t3;
        
        [N_trip_s,C] = hist3([ta(:), tb(:)],'Edges',{ed, ed});
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
        t_from_onset_shuffle{i} = t_from_onset_tmp(:)'; %#ok<AGROW>
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

% triplet vs latency
lat_diff_c = cell(1, n_trial);
lat_peak_diff_c = cell(1, n_trial);
for jj = 1:n_trial
    trip_sig = abs(trip_peak_c{jj}(:,1)) < 100 & abs(trip_peak_c{jj}(:,2)) < 100;
    lat_diff = zeros(n_trip,2);
    for i = 1:n_trip
        lat = latency_c{jj}(ind_trip_c{jj}(i, :));
        lat_a = lat(1) - lat(3);
        lat_b = lat(2) - lat(3);
        lat_diff(i,:) = [lat_a lat_b];
    end
    lat_diff_c{jj} = [lat_diff(trip_sig, 1) lat_diff(trip_sig, 2)];
    lat_peak_diff_c{jj} = [ trip_peak_c{jj}(trip_sig, 1) trip_peak_c{jj}(trip_sig, 2)];

end


% 
% a = cell2mat(lat_peak_diff_c(:));
% b = cell2mat(lat_diff_c(:));
% figure(5)
% plot(a(:),b(:),'.')

R.triplet.sh_ind_c = sh_ind;
R.triplet.PETH_1st_c = PETH_1st_c;
R.triplet.PETH_2nd_c = PETH_2nd_c;
R.triplet.latency_1st_c = latency_1st_c;
R.triplet.latency_2nd_c = latency_2nd_c;
R.triplet.PETH_c = PETH_c;
R.triplet.latency_c = latency_c;
R.triplet.PETH_unique_c = PETH_unique_c;
R.triplet.lat_diff_c  = lat_diff_c;
R.triplet.lat_peak_diff_c = lat_peak_diff_c;
R.triplet.ind_trip_c = ind_trip_c;
R.triplet.N_trip_c = N_trip_c;
R.triplet.N_trip_bin = C;
R.triplet.trip_peak_c = trip_peak_c;
R.triplet.N_trip_shuffle_c = N_trip_shuffle_c;
R.triplet.trip_peak_shuffle_c = trip_peak_shuffle_c;
R.triplet.t_from_onset_shuffle_c =  t_from_onset_shuffle_c;
R.triplet.t_from_onset_c =  t_from_onset_c;
R.triplet.up_onset = up_onset_c;
R.triplet.up_offset = up_offset_c;
disp('Done.');
end
