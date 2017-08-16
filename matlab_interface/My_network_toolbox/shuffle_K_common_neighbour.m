function [ K ] = shuffle_K_common_neighbour( K, I, J, cn_scale_weight )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
show_wait_bar = 0;

fprintf('shuffle_K_common_neighbour....')
if cn_scale_weight > 1
    
    
    %%%%% find number of common pre-synaptic neighbours
    A = sparse(I,J, ones(size(I))); % I:pre, J:post
    cn = A'* A; % A * A' gives common post-synaptic neighbours, note that cn is symmetric
    cn(logical(eye(size(cn)))) = 0; % no self-connection allowed
    
    N = max(max(I), max(J));
    
    K_minmax = minmax(K(:)');
    % now shuffle the K according to common neighbour rule
    cn_minmax = minmax(cn(:)');
    
    if show_wait_bar == 1
        wb_h = waitbar(0,'Please wait...'); wb_p = 0; wb_p_tot = N;
    end
    
    for j = 1:N
        if mod(j*10,round(N/10)*10) == 0
            % fprintf('%d...', 10 - j*10 / (round(N/10)*10));
        end
        if show_wait_bar == 1
            wb_p = wb_p + 1; waitbar(wb_p / wb_p_tot, wb_h);
        end
        c_tmp = cn(I(J == j), j);
        c_tmp = c_tmp(:)';
        K_tmp = K(J == j);
        K_tmp = K_tmp(:)';
        K_tmp_shuffle = zeros(size(K_tmp));
        % bilinear factor: the following scheme not doing much because ck_factor is
        % around 1.1, a very small value due to using global cn_minmax and
        % K_minmax
        ck_factor = 1 + (c_tmp - cn_minmax(1))./diff(cn_minmax) .* ...
            (K_tmp - K_minmax(1))./diff(K_minmax) * (cn_scale_weight - 1);
        max(ck_factor)
        [~, ck_factor_sort_ind] = sort( rand(size(c_tmp))./ck_factor, 'ascend' );
        % shuffle
        K_tmp_shuffle(ck_factor_sort_ind) = sort( K_tmp, 'descend' );
        K(J == j) =  K_tmp_shuffle;
    end
    if show_wait_bar == 1
        close(wb_h)
    end
    fprintf('\n');
    
elseif cn_scale_weight == 1
    disp('cn_scale_weight == 1, shuffling skipped.')
elseif cn_scale_weight < 1
    disp('cn_scale_weight cannot be less that 1.')
end
end

