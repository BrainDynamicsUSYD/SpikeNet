
%%%%%%%%%%%% load data here
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
% The above table uses tau_m = 10 instead of 15!!!!!!!
% Generate a new table???
% Is tau_m's effect on the table linear??

% synI{}.I_mean
% popV.V_mean{}
V_mean = R.popV.V_mean{1};
V_std = R.popV.V_std{1};
I_mean = R.synI{1}.I_mean;
I_std = R.synI{1}.I_std;

step_tot = R.step_tot;

%%%%%%%%%%% compare I_mean as output?
I_mean_s = zeros(1, step_tot); % s for simplified model
I_std_s = zeros(1, step_tot);
for t = 2:step_tot
    mu_V_pre = V_mean(t-1);
    sigma_V_pre =  sqrt(V_std(t-1));
    [ out ] = get_eq3_4_table_lookup(mu_V_pre, sigma_V_pre, table,  'both');
    I_mean_s(t) = out(1); % proportional to
    I_std_s(t) = out(1)*out(2)^2; % proportional to
end


% compare V_mean time serious?









% solve the mean and CV of first passage time of the O-U process, i.e, (eq 3.4)
