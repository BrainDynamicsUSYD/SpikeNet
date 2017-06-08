% R = load('0001-201705081807-23570_in_1494231073516_out_RYG.mat','reduced');
R =load('0002-201705081807-23676_in_1494231075090_out_RYG.mat','reduced');

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
n_trial = 36;
n_circle = 1000;
n_sample = 60;
L = lattice_nD(2, 31);
% I_sample = find(  L(:,1).^2+L(:,2).^2 <= n_circle/pi ); % circular method
I_sample = find(  L(:,1).^2 <= 15^2 ); % rectangular method
sh0 = sh0(I_sample,:);


% for jj = 1:n_trial % need this to get a wide range of latency

% sampling neurons
sh = sh0(randperm(length(I_sample),n_sample), :);

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

figure(1)
hold on;
plot(latency_1st, latency_2nd,'.')
[C,P] = corrcoef(latency_1st,latency_2nd,'rows','complete') %#ok<*NOPTS>

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

mean(PETH_unique)
median(PETH_unique)
std(PETH_unique)

% figure(2)
% for i = 1:36
%     subplot(6,6,i); hold on;
%     plot(PETH_1st(i,:),'r');
%     plot(PETH_2nd(i,:),'b');
% end

% triplet
figure(3)
trip_peak = zeros(n_trial, 3);
ind_trip_mat = zeros(n_trial, 3);
for i = 1:n_trial
    i 
    ah(i) = subaxis(6,6,i,'PR',0.01);
    ind_trip = randperm(n_sample, 3);
    
    ind_trip_mat(i, :) = ind_trip(:)';
    sp = sh(ind_trip, :);
 
    [t1, t2, t3] = meshgrid( find(sp(1,:)), find(sp(2,:)), find(sp(3,:)) );
    ta = t1-t3;
    tb = t2-t3;
    
   
    
    [N,C] = hist3([ta(:), tb(:)],'Edges',{ed, ed});
    N = N / (sum(sp(1,:)) * (bin/1000)^2 );
    N = imgaussfilt(N,gauss_sigma_bin);
    imagesc(C{1},C{2}, N)
    axis([-150 150 -150 150])
    
    [n_max,n_ind] = max(N(:));
    [n_i,n_j] = ind2sub(size(N),n_ind);
    trip_peak(i,:) = [C{1}(n_i), C{2}(n_j), n_max];

end
set(gca,'CLIM',[0 320])
colormap('jet')

% shuffled triplet 
trip_peak_shuffle = zeros(n_trial, 3);
sh_shuffle = raster_marginals_shuffling(sh);
figure(4);
for i = 1:n_trial
    i 
    ah(i) = subaxis(6,6,i,'PR',0.01);
    ind_trip = ind_trip_mat(i,:);
    
    sp = sh_shuffle(ind_trip, :);

    [t1, t2, t3] = meshgrid( find(sp(1,:)), find(sp(2,:)), find(sp(3,:)) );
    ta = t1-t3;
    tb = t2-t3;
    
    [N,C] = hist3([ta(:), tb(:)],'Edges',{ed, ed});
    N = N / (sum(sp(1,:)) * (bin/1000)^2 );
    N = imgaussfilt(N,gauss_sigma_bin) ;
    imagesc(C{1},C{2}, N)
    axis([-150 150 -150 150])
    
    [n_max,n_ind] = max(N(:));
    [n_i,n_j] = ind2sub(size(N),n_ind);
    trip_peak_shuffle(i,:) = [n_i,n_j,n_max];
end

set(gca,'CLIM',[0 320])
colormap('jet')

% triplet vs latency
trip_sig = (trip_peak(:,3)./trip_peak_shuffle(:,3) > 2);
lat_diff = zeros(n_trial,2);
for i = 1:n_trial
    lat = latency_1st(ind_trip_mat(i, :));
    lat_a = lat(1) - lat(3);
    lat_b = lat(2) - lat(3);
    lat_diff(i,:) = [lat_a lat_b];
end

figure(5); hold on;
plot(lat_diff(trip_sig, 1), trip_peak( trip_sig, 1), '.' )
plot(lat_diff(trip_sig, 2), trip_peak(trip_sig, 2), '.' )




