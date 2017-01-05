function [ x, y, width, mlh ] = fit_bayesian_bump_2_spikes( spike_x_pos_o,spike_y_pos_o, fw, mode)
%UNTITLED2 Summary of this function goes her
%   Detailed explanation goes here

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
    
    
    n_std = 4;
    x_r = ((-n_std*(ceil(width_guess))):(n_std*(ceil(width_guess)))) + round(mean(spike_x_pos_shifted));
    y_r = ((-n_std*(ceil(width_guess))):(n_std*(ceil(width_guess)))) + round(mean(spike_y_pos_shifted));
    [x_grid, y_grid] = meshgrid( x_r, y_r);
    spike_count = hist3([spike_x_pos_shifted spike_y_pos_shifted],'Edges',{x_r,y_r});
    peak_rate_guess = 2*max(max(spike_count));
    
    x0 = [x_c_guess, y_c_guess, peak_rate_guess, width_guess ];
    
    func_handle = @(v) -get_lh_log_ad_hoc(x_grid, y_grid, spike_count, v(1), v(2), v(3), v(4) );
    
    switch lower(mode)
        case 'bayesian'
            [v, fval] = fminsearch(func_handle, x0);
            mlh = -fval;
        case 'quick'
            v = x0;
            mlh = 0;
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


end

function lh_log = get_lh_log_ad_hoc(x_grid, y_grid, spike_count, x_c, y_c, peak_rate, width )

peak_rate = abs(peak_rate); % should always be positive
gauss_rate =  exp( -((x_grid - x_c).^2 + (y_grid - y_c).^2)/(2*width^2) );
gauss_rate = gauss_rate/max(max(gauss_rate))*peak_rate;
lh_log = sum(sum( spike_count.*log(gauss_rate) - gauss_rate ));

end
