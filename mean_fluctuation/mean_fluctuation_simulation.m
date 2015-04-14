function mean_fluctuation_simulation(varargin)


% varargin is for PBS arrary job
if nargin == 0
    clc;clear all;close all;
    cd /import/yossarian1/yifan/Project1/
    addpath(genpath(cd));
    cd tmp_data
end

loop_num = 0;

for loop = 1:50;
    
    loop_num = loop_num + 1;
    
    % For PBS array job
    if nargin ~= 0
        PBS_ARRAYID = varargin{1};
        if loop_num ~=  PBS_ARRAYID
            continue;
        end
    end
    
    % seed the matlab rand function! The seed is global.
    % Be very careful about that!!!!!!!!!!!
    date_now = datestr(now,'yyyymmddHHMM-SSFFF');
    scan_temp = textscan(date_now,'%s','Delimiter','-');
    rand_seed = loop_num*10^5+eval(scan_temp{1}{2}); % scan_temp{1}{2} is the SSFFF part
    rng(rand_seed,'twister'); %  Note that this effect is global!
    fprintf('Random number generator seed is: %f\n', rand_seed );
    
    
    
    % time
    dt = 0.1; % ms
    step_tot = 100*10^4; % 100 bio-second took 70mins to simulate
    % T = (1:step_tot)*dt/1000;% sec
    
    
    % %% define the parameters
    t_m = 10; %ms
    poss = 1;
    for c_mu = 0.1; %5
        for  c_sigma = 21 % 20.5
            if ~isreal(4*c_mu^2 - 8*(c_mu^2-c_sigma^2))
                disp('Following combination is not possible!');
                c_mu, c_sigma
                poss = 0;
                
            end
        end
    end
    
    
    
    if poss == 1
        %% load pre-calculated table for (eq 3.4)
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
        
        
        
        %% simulation
        
        N_pop = 2;
        % c is number of connection
        c = [1 1;    % [EE  IE
            1 1];   %  EI  II]
        % j is the connection strength
        j_E = (2*c_mu + sqrt(4*c_mu^2 - 8*(c_mu^2-c_sigma^2)))/4;
        j_I = j_E - c_mu;
        
        
        j = [j_E j_E;    % [EE  IE
            -j_I -j_I];   %  EI  II]
        % j_EI can be changed to -13.56 and the results are still similar
        
        % the above c and j values give c_miu = 5 and c_sigma = 20.2
        
        % bookkeeping
        miu_V = zeros(N_pop, step_tot);
        sigma2_V = zeros(N_pop, step_tot);
        firing_rate = zeros(N_pop, step_tot);
        CV_hist =  zeros(N_pop, step_tot);
        % initial condition
        miu_V(:,1) = ones(N_pop,1)*5.124;
        sigma2_V(:,1) = ones(N_pop,1)*5.496^2;
        
        
        % external inputs
        %%%%%%%% try two correlated noises for mu_ext and sigma_ext

        tau_1 = 100; % what's it in the full system?
        r1 = exp_corr_gaussian_noise(step_tot,tau_1);
        
        tau_2 = 100; %
        r_mid = exp_corr_gaussian_noise(step_tot,tau_2);
        
        alpha = 0.4; % zero time-lag correlation of the two signals
        % what's it in the full system?
        r2 = alpha*r1 + sqrt(1-alpha^2)*r_mid;
        r2 = (r2 - mean(r2))/std(r2);
        % corrcoef(r1, r2)
        
        
        mu_ext = 4.5 + r1; %
        sigma_ext = 5.0 + r2; %
        %%%%%%%% try two correlated noises for mu_ext and sigma_ext
        
        
        tic;  % 4e4 steps take about 1 min to finish
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
            
        end
        toc;
        
    end
    
    
    filename = sprintf('%04d-mean-fluc', loop_num);
    save(filename, 'miu_V', 'sigma2_V', 'firing_rate', 'CV_hist','mu_ext', 'sigma_ext',...
        'dt', 'step_tot', 'c_mu','c_sigma', 't_m', 'tau_1','tau_2', 'alpha');
    
end


end




