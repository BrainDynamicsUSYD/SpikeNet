function [ x, y, width, mlh, height, bayes_factor  ] = fit_bayesian_bump_2_spikes( spike_x_pos_o,spike_y_pos_o, fw, mode)
%UNTITLED2 Summary of this function goes her
%   Detailed explanation goes here
% mode: 'bayesian' or 'quick'


if isempty(spike_x_pos_o)
    x = NaN;
    y = NaN;
    width = NaN;
    mlh = NaN;
    height = NaN;
else
    
    hw = fw/2;
    dd = pi*10^-1;
    spike_x_pos_o = spike_x_pos_o(:);
    spike_y_pos_o = spike_y_pos_o(:);
    
    
    % periodic boundary condition
    x_s = linspace(-hw, hw, 7) + dd;
    y_s = linspace(-hw, hw, 7) + dd;
    [x_s_grid, y_s_grid] = meshgrid(x_s, y_s);
    
    x_shift = x_s_grid(:);
    y_shift = y_s_grid(:);
    width_guess = inf;
    for i = 1:length(x_shift)
        
        % shift
        spike_x_pos_shifted =  mod(spike_x_pos_o - x_shift(i)+hw, fw) - hw;
        spike_y_pos_shifted =  mod(spike_y_pos_o - y_shift(i)+hw, fw) - hw;
        
        width_guess_tmp = max( std(spike_x_pos_shifted), std(spike_y_pos_shifted) );
        
        if width_guess_tmp < width_guess
            width_guess =  width_guess_tmp;
            i_pbc= i;
        end
    end
    
    
    for i = i_pbc
        % shift
        spike_x_pos_shifted =  mod(spike_x_pos_o - x_shift(i)+hw, fw) - hw;
        spike_y_pos_shifted =  mod(spike_y_pos_o - y_shift(i)+hw, fw) - hw;
        x_c_guess = mean(spike_x_pos_shifted);
        y_c_guess = mean(spike_y_pos_shifted);
        width_guess = max( std(spike_x_pos_shifted), std(spike_y_pos_shifted) );
        
        x_r = (-hw:hw) + round(mean(spike_x_pos_shifted));
        y_r = (-hw:hw) + round(mean(spike_y_pos_shifted));
        [x_grid, y_grid] = meshgrid( x_r, y_r);
        spike_count = hist3([spike_x_pos_shifted spike_y_pos_shifted],'Edges',{x_r,y_r});
        peak_rate_guess = 2*max(max(spike_count));
        
        x0 = [x_c_guess, y_c_guess, peak_rate_guess, width_guess ];

        switch lower(mode)
            case 'bayesian'
                func_handle_gauss = @(v) -get_gaus_lh_log_ad_hoc(x_grid, y_grid, spike_count, v(1), v(2), v(3), v(4) );
                [v, fval] = fminsearchbnd(func_handle_gauss, x0, [min(x_r) min(y_r) 0 0] , [max(x_r) max(y_r)  500 100*fw]);
                mlh = -fval;
                
                % get bayes factor
                func_handle_unif = @(rate) -get_unif_lh_log_ad_hoc( spike_count, rate );
                [rate_unif, fval] = fminsearchbnd(func_handle_unif, peak_rate_guess/2, 0 , peak_rate_guess);
                mlh_unif = -fval;
                bayes_factor = mlh - mlh_unif - 0.5*(4-1)*log(sum(sum(spike_count)));
                % See Gu and Gong, 2016, The dynamics of memory retrieval in hierarchical networks: a modeling study

                % get P-value
                % % No!! It is very hard to design a meaningful p-value
                % here.
                
            case 'quick'
                v = x0;
                mlh = 0;
                bayes_factor = NaN;
        end
        
        %     % shift the position back to normal
        %     xx_o = unique(spike_x_pos_o);
        %     xx_shifted = mod(xx_o - x_shift(i)+hw, fw) - hw;
        %     [xx_shifted_sorted, ind_tmp] = sort(xx_shifted);
        %     v(1) = interp1(xx_shifted_sorted,  xx_o(ind_tmp), v(1));
        %
        %     yy_o = unique(spike_y_pos_o);
        %     yy_shifted = mod(yy_o - y_shift(i)+hw, fw) - hw;
        %     [yy_shifted_sorted, ind_tmp] = sort(yy_shifted);
        %     v(2) = interp1(yy_shifted_sorted,  yy_o(ind_tmp), v(2));
        v(1) = v(1) + x_shift(i);
        if v(1) < -hw
            v(1) = v(1) + fw;
        elseif v(1) > hw
            v(1) = v(1) - fw;
        end
        v(2) = v(2) + y_shift(i);
        if v(2) < -hw
            v(2) = v(2) + fw;
        elseif v(2) > hw
            v(2) = v(2) - fw;
        end
        v_final = v;
    end
    
    x = v_final(1);
    y = v_final(2);
    width = v_final(4);
    height = v_final(3);
    
end

end

function lh_log = get_gaus_lh_log_ad_hoc(x_grid, y_grid, spike_count, x_c, y_c, peak_rate, width )

peak_rate = abs(peak_rate); % should always be positive
gauss_rate =  exp( -((x_grid - x_c).^2 + (y_grid - y_c).^2)/(2*width^2) );
gauss_rate = gauss_rate/max(max(gauss_rate))*peak_rate;
lh_log = sum(sum( spike_count.*log(gauss_rate) - gauss_rate ));

end

function lh_log = get_unif_lh_log_ad_hoc( spike_count, rate )
lh_log = sum(sum( spike_count.*log(rate) - rate ));
end
