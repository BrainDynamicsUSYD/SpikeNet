
loop_num = 1;

% seed the matlab rand function! The seed is global.
% Be very careful about that!!!!!!!!!!!
date_now = datestr(now,'yyyymmddHHMM-SSFFF');
scan_temp = textscan(date_now,'%s','Delimiter','-');
rand_seed = loop_num*10^5+eval(scan_temp{1}{2}); % scan_temp{1}{2} is the SSFFF part
rng(rand_seed,'twister'); %  Note that this effect is global!
fprintf('Random number generator seed is: %f\n', rand_seed );

% time
dt = 0.1; % ms
step_tot = 1*10^3; % 100 bio-second took 70mins to simulate
% T = (1:step_tot)*dt/1000;% sec

% %% define the parameters
t_m = 10; %ms    <--- scale it carefully to 15ms !?
poss = 1;
for c_mu = 0.1; %5
    for  c_sigma = 21 % 20.5
        if ~isreal(4*c_mu^2 - 8*(c_mu^2-c_sigma^2))
            disp('The following combination is not possible!');
            c_mu, c_sigma
            poss = 0;
        end
    end
end

% simulation
N_pop = 9;
E_factor = 1.5;

% c is number of connection
%         c = [1 1;    % [EE  IE
%             1 1];   %  EI  II]

c = ones(N_pop, N_pop);

% j is the connection strength
j_E = (2*c_mu + sqrt(4*c_mu^2 - 8*(c_mu^2-c_sigma^2)))/4;
j_I = j_E - c_mu;

%         j = [j_E j_E;    % [EE  IE
%             -j_I -j_I];   %  EI  II]

lesion_1 = 1.1;
lesion_2 = 1;
lesion_3 = 1;
lesion_4 = 0.6;
[P, CL] = inter_module_Pmatrix(500*ones(8, 1), 0.2, 0.7);
P(CL==1) = P(CL==1)*lesion_1;
P(CL==2) = P(CL==2)*lesion_2;
P(CL==3) = P(CL==3)*lesion_3;
P(CL==4) = P(CL==4)*lesion_4;
P = P/sum(P(1,:))*E_factor;

j = j_E*ones(N_pop, N_pop)/8;
j(1:8,1:8) = P*j_E;
j(end, :) = -j_I;

% j_EI can be changed to -13.56 and the results are still similar

% the above c and j values give c_miu = 5 and c_sigma = 20.2

% simulate_mean_fluc

mu_ext = 2*ones(1, step_tot);

[ R ] = simulate_mean_fluc(N_pop, step_tot, dt, c, j, mu_ext );

figure(loop_num);
for i = 1:N_pop
    subplot(N_pop, 1, i);
    plot(dt*(1:step_tot), R.firing_rate(i,:));
end







