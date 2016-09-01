
function visualize_grid_and_SWR(R, LFP_centre_x, LFP_centre_y)
close all;
clc;

ind_a_vec = R.grid.raw.ind_ab(1,:);
ind_b_vec = R.grid.raw.ind_ab(2,:);

hw = 31;
fw = 2*hw+1;
x_mean_chosen = R.grid.centre(1,:);
y_mean_chosen = R.grid.centre(2,:);
dist_std_chosen = R.grid.radius;
t_mid = R.grid.raw.t_mid;
t_mid_chosen = R.grid.t_mid;
jump_dist = R.grid.jump_dist;

[Lattice, ~] = lattice_nD(2, hw);
x_pos_o = Lattice(R.spike_hist_compressed{1}, 1);
y_pos_o = Lattice(R.spike_hist_compressed{1}, 2);


% %
% LFP_centre_x = linspace(-hw, hw, 5);
% LFP_centre_y = linspace(-hw, hw, 5);
% [LFP_centre_x, LFP_centre_y] = meshgrid(LFP_centre_x(1:end-1), LFP_centre_y(1:end-1));
% LFP_centre_x = LFP_centre_x(:);
% LFP_centre_y = LFP_centre_y(:);

figure('Name','Vis','color','w','NumberTitle','off');

axis equal;
box on;
set(gca,'xtick',[],'ytick',[]);

xlim([-hw hw]);
ylim([-hw hw]);
hold on;


%
plot(LFP_centre_x, LFP_centre_y, '.b', 'MarkerSize',69);

j = 1;

ang=0:0.01:2*pi;
x_shift_vs = [0 fw fw -fw -fw fw -fw 0 0 ];
y_shift_vs = [0 fw -fw fw -fw 0  0   fw -fw];
win_len = R.grid.raw.win_len;
for t = 1:length(t_mid);
    h2 = plot(100,0);
    h3 = plot(100,0);
    

    t_range_tmp = (t_mid(t)-win_len/2):(t_mid(t) + win_len/2 -1);
    % ripple dectection
    detceded = sum(R.LFP.ripple_event.is_SWR(:, t_range_tmp), 2) > 0;
    h4 = plot(LFP_centre_x(detceded), LFP_centre_y(detceded), '.r', 'MarkerSize',69);
    
    ind_range_tmp = ind_a_vec(t):ind_b_vec(t);
    h1 = plot(x_pos_o(ind_range_tmp), y_pos_o(ind_range_tmp), 'bo');
    
    if sum(t_mid_chosen == t_mid(t)) == 1
        
        x_tmp = x_mean_chosen(j);
        y_tmp = y_mean_chosen(j);
        r_cos = x_tmp+dist_std_chosen(j)*cos(ang);
        r_sin = y_tmp+dist_std_chosen(j)*sin(ang);
        
        h3 = plot(r_cos - x_shift_vs(1),r_sin - y_shift_vs(1),'b', ...
            r_cos - x_shift_vs(2),r_sin - y_shift_vs(2),'b',...
            r_cos - x_shift_vs(3),r_sin - y_shift_vs(3),'b',...
            r_cos - x_shift_vs(4),r_sin - y_shift_vs(4),'b',...
            r_cos - x_shift_vs(5),r_sin - y_shift_vs(5),'b', ...
            r_cos - x_shift_vs(6),r_sin - y_shift_vs(6),'b', ...
            r_cos - x_shift_vs(7),r_sin - y_shift_vs(7),'b', ...
            r_cos - x_shift_vs(8),r_sin - y_shift_vs(8),'b', ...
            r_cos - x_shift_vs(9),r_sin - y_shift_vs(9),'b');
        
        if j == 1
            h2 = plot( x_tmp, y_tmp, 'r>', 'MarkerSize', 8);
        else
            h2 = plot( x_tmp, y_tmp, 'rx', 'MarkerSize', 16);
            xs = sprintf('jump_size = %4.2f, radius = %4.2f', jump_dist(j-1), dist_std_chosen(j));
            xlabel(xs);
        end
        
        
        j = j + 1;
    end
    pause(0.1);
    
    
    delete(h1);
    delete(h2);
    delete(h3);
    delete(h4);
    ts = sprintf('time = %8.1f ms', t_mid(t)*0.1);
    title(ts);
    
    
end
end
