function [ K ] = get_K_given_certain_conditions( I, J, K_mu_log, K_sigma_log, K_neuron_sum, plot_results)
%[ K ] = get_K_given_certain_conditions( I, J, K_mu_log, K_sigma_log, K_neuron_sum, plot_results)
%   This function returns one stochastic soluton of the following problem:
%   Given (1) the distribution of synaptic weights pooled from all of the
%   neurons (2) the sum of incoming synaptic weights for each neuron and
%   (3) the in-degree of each neruon, find the incoming synaptic weights
%   for each neuron that satisfy all of the above 3 conditions.
% 
%   This problem is non-trivial.
%   My solution is in publishable form. But is it publishable?
%
%   Yifan Gu, School of Physics, USYD, Mar 2016.
%
%
% %%%% example code for using this function:
%
% hw = 11;
% P0_init = 0.1;
% degree_CV = 0.5;
% tau_c = 10;
% in_out_corrcoef = 0.2;
% 
% [ll_matrix, N_e, ~] = spatial_embed_in_out_network(hw, P0_init, degree_CV, tau_c, in_out_corrcoef);
% [I_e,J_e,~] = find(ll_matrix);
% A = sparse(I_e,J_e, ones(size(I_e)));
% clear ll_matrix;
% 
% K11 = 0.5e-3; % miuSiemens
% EE_CV = 2;
% K_mu_log = K11; %/(mean(in_degree));
% K_sigma_log = K_mu_log*EE_CV;
% 
% % define and in_degree vs K_neuron_sum scaling 
% K_tot = K_mu_log*length(I_e);
% in_degree = full(sum(A, 1));
% scale = sqrt(in_degree);
% K_neuron_sum = scale*K_tot/sum(scale);
% 
% plot_results = 1;
% tic
% [ K ] = get_K_given_certain_conditions( I_e, J_e, K_mu_log, K_sigma_log, K_neuron_sum, plot_results);
% toc
%
%%%% example code for using this function ends here


in_degree = full(sum(sparse(I,J, ones(size(I))), 1));

N = max( length(unique(I)), length(unique(J)) ); % number of neurons
K_length = length(I);
K = zeros(1, K_length);

Mu_log = K_mu_log;
Sigma_log = K_sigma_log;
Mu_norm = log((Mu_log.^2)./sqrt(Sigma_log.^2+Mu_log.^2));
Sigma_norm = sqrt(log(Sigma_log.^2./(Mu_log.^2)+1));

% K_sum = Mu_log*K_length;
% scale = sqrt(in_degree);
% K_neuron_sum = scale*K_tot/sum(scale);

% check input condition consistency
K_tot = sum(K_neuron_sum);
if abs(K_tot - Mu_log*K_length)/K_tot > 0.05
    warning('The given K_neuron_sum is not consistent with K_mu_log!')
end

err_max = 0.02; % unit: percentage
max_while_count = 100;

j_left = 1:N;
std_K_neuron = zeros(1,N);
K_pool_num = 0; % usually needs a 2nd pool
while sum(isnan(j_left)) < N
    K_pool_num = K_pool_num + 1;
    % generate the pool for K
    K_pool_left = exp(randn(1, K_length)*Sigma_norm + Mu_norm);
    K_pool_left = sort(K_pool_left); % from small values to big values
    
    err_acc = 0;
    % find suitable weights for each neuron from the current pool if
    % possible
    for jj = j_left
        if ~isnan(jj) % if not found in the previous pool
            err = 1;
            err_sign_match = true;
            while_count = 0;
            while while_count <  max_while_count && ...
                    err > err_max || err_sign_match
                
                N_s = -1; % s for smaller
                N_b = -1; % b for bigger
                N_s_tot = -2;
                N_b_tot = -2;
                
                % find a random point that separates the current pool into a
                % smaller and a larger one
                while while_count <  max_while_count && ...
                        (N_s < 0 || N_b < 0 || N_s > N_s_tot || N_b > N_b_tot) % loop until the point meets the criteria
                    
                    bs_sep = randperm(length(K_pool_left)-1,1);
                    % find the right number of samples from the two
                    % sub-pools
                    A = [mean(K_pool_left(1:bs_sep)) mean(K_pool_left(bs_sep+1:end));
                        1  1];
                    B = [K_neuron_sum(jj); in_degree(jj)];
                    N_sb = linsolve(A,B); %X = linsolve(A,B) solves the linear system A*X=B
                    N_s = round(N_sb(1));
                    N_b = in_degree(jj) - N_s;
                    N_s_tot = bs_sep;
                    N_b_tot = length(K_pool_left) - bs_sep;
                    while_count = while_count + 1;
                end
                % Control the maximum number of attempts at finding the
                % suitable separation point.
                % If the max_while_count is exceeded, try it in the
                % next sample pool
                if while_count == max_while_count
                    break;
                end
                % sample from the two sub-pools
                sub_s = randperm(N_s_tot, N_s);
                sub_b = randperm(N_b_tot, N_b) + N_s_tot;
                
                % control the accummulated error (err_acc)
                % when the previous accummulated error is positive, prefer
                % negative error time round, and vice versa, so that it
                % stays around zero.
                K_neuron_tmp = K_pool_left([sub_s sub_b]);
                err = abs( (sum(K_neuron_tmp) - K_neuron_sum(jj))/K_neuron_sum(jj) );
                err_sign = sign(sum(K_neuron_tmp) - K_neuron_sum(jj));
                err_sign_match = ( err_sign == sign(err_acc) );
            end
            
            if while_count < max_while_count
                j_left(jj) = NaN; % NaN means the weights for neuron #j have been successully sampled
                K(J == jj) = K_pool_left([sub_s sub_b]);
                std_K_neuron(jj) = std( K_pool_left([sub_s sub_b])); % record std, which is not pre-defined nor controlled
                K_pool_left([sub_s sub_b]) = [];
                err_acc= err_acc + sum(K_neuron_tmp) - K(jj);
                % clc; P=P+1; disp(P); % show progress
            end
            
        end
        
    end
end

%%% show results
if nargin < 6
    plot_results = 0;
end
if plot_results == 1
    figure('NumberTitle','off','name','get_K_given_in_degree','color','w')
    % assume log-normal K distribution
    subplot(2,1,1);
    EPSP = g_2_EPSP( K ); % log-normal K distribution leads to log-normal EPSP distribution??
    [Y, X] = hist(log10(EPSP), 100);
    f = fit(X(:), Y(:),'gauss1');
    plot(f,X,Y)
    ylabel('count')
    xlabel(sprintf('log10(EPSP) (mV), mean = %f, max = %f.', mean(EPSP),max(EPSP)))
    % compare generated K and desired K (sum for each neuron)
    K_neuron_sum_gen = sum(sparse(I, J, K),1); % gen for generated
    subplot(2,1,2);hold on;box on;
    errorbar(in_degree, K_neuron_sum_gen, std_K_neuron, 'bo');
    [in_degree_s, ind_s] = sort(in_degree);
    plot(in_degree_s, K_neuron_sum(ind_s), 'r');
    xlabel('in degree')
    ylabel('summed synaptic weights')
    legend('generated','desired')
    % %
    % subplot(3,1,3);hold on;
    % for j = 1:N_e
    %     plot(j*ones(1,sum(J_e == j)), K_e(J_e == j), 'bx');
    % end
end

end
