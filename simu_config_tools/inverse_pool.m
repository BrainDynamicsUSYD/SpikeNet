function [ K ] = inverse_pool(  K_num, K_scale, pool_generator)
% [ K ] = inverse_pool(  K_num, K_scale, pool_generator)
% K_num should be a vector that specifies the number of elements in each sub-pool K{i}
% K_scale should specify the (relative ratio of) sum of the elements in each sub-pool K{i}
% pool_generator should be a function handle that generates the total pool
% The output K is a cell contains the sub-pools K{i} that satisfy the given
% conditions.

% % Example
% Mu_norm = log((Mu_log.^2)./sqrt(Sigma_log.^2+Mu_log.^2));
% Sigma_norm = sqrt(log(Sigma_log.^2./(Mu_log.^2)+1));
% K_pool_generator = @(N)exp(randn(1, N)*Sigma_norm + Mu_norm);


% parameters
show_wait_bar = 0;
max_K_pool_num = 10; % ?
err_max = 0.1; % unit: percentage
max_sb_round = 4; % ?
max_err_round = 50; % ?
pool_length_extra = 0.01; % ??


pool_length = round(sum(K_num)*(1+pool_length_extra));
% K_sum specify the  sum of the elements in each sub-pool K{i}
K_sum_tot = sum(pool_generator(sum(K_num)));
K_sum = K_scale./sum(K_scale)*K_sum_tot;
N = length(K_num); % number of slots

% re-organize j_left, K_num and K_sum to prioritize small and
% large K_num requirements
[~, K_num_ind] = sort(K_num);
N_h = round(N/2);
AZBY = zeros(1,N);
AZBY(1:2:end) = 1:N_h;
AZBY(2:2:end) = N:-1:N_h+1;
K_num_AZBY = K_num( K_num_ind( AZBY ) );
i_left = 1:N;
i_left_AZBY = i_left(K_num_ind( AZBY ));
K_sum_AZBY = K_sum(K_num_ind( AZBY ));


K = cell(1,N);
K_std = zeros(1,N);
K_pool_num = 0;
K_pool_left = [];
err_acc = 0;
while sum(isnan(i_left_AZBY)) < N &&  K_pool_num < max_K_pool_num


    K_pool_num = K_pool_num + 1;
    fprintf('Generate K: K-pool number = %d.\n ', K_pool_num)
    
    % generate the pool for K, recycle the previous pool_left
    K_pool = pool_generator(pool_length);
    K_pool_left = [K_pool_left K_pool ] ; %#ok<AGROW>
    K_pool_left = sort(K_pool_left); % from small values to big values

    % quick check for solvability
    solvable = true(1,N);
    for i_check = 1:N
        if sum(K_pool_left(1:K_num(i_check))) >  K_sum(i_check) ...
                || sum(K_pool_left(end-K_num(i_check)+1:end)) <  K_sum(i_check)
            solvable(i_check) = false;
        end
    end
    N_unsolvable = sum( ~solvable );
    
    if  N_unsolvable > 0 &&  N_unsolvable < 0.05*N
        fprintf('There are %d (K_num, K_sum) pair cannot be solved.\n', N_unsolvable);
        fprintf('Trying another pool...\n');
        K_pool_left = [];
        continue;
    elseif N_unsolvable >= 0.05*N
        fprintf('There are %d (K_num, K_sum) pair cannot be solved.\n', N_unsolvable);
        fprintf('Inverse-pooling process abandoned. Adjust the input arguments!\n');
        break;
    else
        % find suitable weights for each neuron from the current pool if
        % possible
        if show_wait_bar == 1
            wb_h = waitbar(0,'Please wait...');
            wb_p = 0;
            wb_p_tot = nnz(i_left_AZBY);
        end
        for ind = 1:N
            ii = i_left_AZBY(ind);
            
            if show_wait_bar == 1
                wb_p = wb_p + 1;
                waitbar(wb_p / wb_p_tot, wb_h);
            end
            
            if ~isnan(ii) % if not found in the previous pool
                sb_round = 0;
                K_sum_tmp = K_sum_AZBY(ind);
                K_num_tmp = K_num_AZBY(ind);
                while  sb_round < max_sb_round
                    sb_round = sb_round + 1;
                    % find suitable separation point
                    [ N_s_tot, N_b_tot, N_s, N_b ] = find_bs_sep(K_pool_left, K_sum_tmp, K_num_tmp);
                    if isnan(N_s_tot) % if cannot find, try it in the next K-pool
                        fprintf('%d is unsolvable.\n', ii);
                        break;
                    else
                        err_round = 0;
                        while err_round <  max_err_round
                            err_round = err_round + 1;
                            % sample from the two sub-pools
                            sub_s = randperm(N_s_tot, N_s);
                            sub_b = randperm(N_b_tot, N_b) + N_s_tot;
                            
                            % control the accummulated error (err_acc)
                            % when the previous accummulated error is positive, prefer
                            % negative error time round, and vice versa, so that it
                            % stays around zero.
                            K_selected_tmp = K_pool_left([sub_s sub_b]);
                            err = abs( (sum(K_selected_tmp) - K_sum_tmp)/K_sum_tmp );
                            err_sign = sign(sum(K_selected_tmp) - K_sum_tmp);
                            err_sign_match = ( err_sign == sign(err_acc) );
                            
                            % if the error is good, stop the while loop
                            if err <= err_max && ~err_sign_match
                                i_left_AZBY(ind) = NaN; % NaN means the weights for neuron #j have been successully sampled
                                K{ii} = K_selected_tmp;
                                K_std(ii) = std( K_selected_tmp ); % record std, which is not pre-defined nor controlled
                                K_pool_left([sub_s sub_b]) = [];
                                err_acc= err_acc + sum(K_selected_tmp) - K_sum_tmp;
                                
                                sb_round = max_sb_round; % help to break out of the outer loop!
                                break;
                            end
                        end % error_round
                    end
                end % sb_round
            end
        end
        if show_wait_bar == 1
            close(wb_h)
        end
    end
 
    fprintf('%d out of %d reverse-pooled.\n', sum(isnan(i_left_AZBY)), N);
end


if sum(isnan(i_left_AZBY)) < N
    K = {NaN};
    fprintf('Inverse-pooling process abandoned. Adjust the input arguments!\n');
end


end

function [ N_s_tot, N_b_tot, N_s, N_b ] = find_bs_sep(K_pool_left, K_sum, K_num)
% this function finds the proper separation point
bs_sep_found = false;

while_count = 0;
max_while_count = 2*log2(length(K_pool_left)); %

pool_length = length(K_pool_left);

sep_lower_limit = 1;
sep_upper_limit = length(K_pool_left)-1;
% find a random point that separates the current pool into a
% smaller and a larger one
while while_count <  max_while_count && sep_upper_limit - sep_lower_limit > 0 ...
        && pool_length >= K_num % loop until the proper separation point is found
    
    while_count = while_count + 1;
    adjust_limits = 0;
    bs_sep = randperm(sep_upper_limit - sep_lower_limit + 1, 1) + (sep_lower_limit - 1);
    
    % find the right number of samples from the two
    % sub-pools
    mean_s = mean(K_pool_left(1:bs_sep));
    mean_b = mean(K_pool_left(bs_sep+1:end));

    if mean_s > K_sum/K_num
       sep_upper_limit = bs_sep - 1;
    elseif mean_b < K_sum/K_num
       sep_lower_limit = bs_sep + 1;
    else
        N_s_tot = bs_sep;% s for smaller
        N_b_tot = pool_length - N_s_tot;% b for bigger
        
        A = [mean_s mean_b;
            1  1];
        B = [K_sum; K_num];
        N_sb = linsolve(A,B); %X = linsolve(A,B) solves the linear system A*X=B
        N_s = N_sb(1);
        N_b = N_sb(2);

        % adjust the searching boundary
        if N_s > N_s_tot
            sep_lower_limit = bs_sep + 1;
            adjust_limits = 1;
        end
        if N_b > N_b_tot
            sep_upper_limit = bs_sep - 1;
            adjust_limits =  adjust_limits + 0.2;
        end
        
        % check if the proper separation point is found
        if adjust_limits == 0
            if N_s < N_b
                N_s = round(N_sb(1));
                N_b = K_num - N_s;
            else
                N_b = round(N_sb(2));
                N_s = K_num - N_b;
            end
            bs_sep_found = true;
            break;
        end
    end
    
end


% if the proper separation point is not found, output NaN
if ~bs_sep_found
    N_s_tot = NaN;
    N_b_tot = NaN;
    N_s = NaN;
    N_b = NaN;
end


end

