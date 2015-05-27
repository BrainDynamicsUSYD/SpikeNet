
function  [ R ] = simulate_mean_fluc(N_pop, step_tot, dt, c, j, mu_ext, sigma_ext )

h = waitbar(0,'Please wait...');

% load pre-calculated table for (eq 3.4)
table = load('v_CV_tables.mat');
% manually fix some (possibly) numerical errors
table.CV_mat(1:13,1) = 1;
table.CV_mat(end-2:end,1) = 0;
table.v_mat(end-2:end,1) = 0.5;
% something is still wrong here!!!
% transpose for interp2
table.miu_V_mat = table.miu_V_mat';
table.sigma_V_mat = table.sigma_V_mat';
table.CV_mat = table.CV_mat';
table.v_mat = table.v_mat';

% memberane time constant
t_m = 10; %ms

% bookkeeping
miu_V = zeros(N_pop, step_tot);
sigma2_V = zeros(N_pop, step_tot);
firing_rate = zeros(N_pop, step_tot);
CV_hist =  zeros(N_pop, step_tot);

% initial condition
miu_V(:,1) = rand(N_pop,1)*5.124;
sigma2_V(:,1) = rand(N_pop,1)*5.496^2;


% external inputs (identical to all populations)
% mu_ext = 0.5*ones(1, step_tot);
if nargin == 6
    sigma_ext = zeros(1, step_tot);
end

% 4e4 steps take about 1 min to finish
for t = 2:step_tot
    
    miu_I = ones(N_pop,1)*mu_ext(t)/t_m;
    sigma2_I = ones(N_pop,1)*sigma_ext(t)^2/(t_m/2);
    % calculate the mean and std of currents
    for i_pre = 1:N_pop
        % solve the mean and CV of first passage time of the O-U process, i.e, (eq 3.4)
        miu_V_pre = miu_V(i_pre, t-1);
        sigma_V_pre =  sqrt(sigma2_V(i_pre, t-1));
        [ out ] = get_eq3_4_table_lookup(miu_V_pre, sigma_V_pre, table,  'both');
        v_pre = out(1);
        CV_pre = out(2);
        firing_rate(i_pre, t) = v_pre;
        CV_hist(i_pre, t) = CV_pre;
        for j_post = 1:N_pop
            c_tmp = c(i_pre, j_post);
            j_tmp = j(i_pre, j_post);
            miu_I(j_post) = miu_I(j_post) + c_tmp*j_tmp*v_pre;
            sigma2_I(j_post) = sigma2_I(j_post) + c_tmp*j_tmp^2*v_pre*CV_pre^2;
        end
    end
    % integrate the ODEs
    d_miu_V = - miu_V(:,t-1)/t_m + miu_I;
    d_sigma2_V = - sigma2_V(:,t-1)/(t_m/2) + sigma2_I;
    
    miu_V(:,t) = miu_V(:,t-1) + d_miu_V*dt;
    sigma2_V(:,t) = sigma2_V(:,t-1) + d_sigma2_V*dt;
    
    % show progress
    if mod( t, 100 ) == 0
        waitbar(t/step_tot,h)
    end
    
end


close(h);
delete(h);

R.firing_rate = firing_rate;
end
