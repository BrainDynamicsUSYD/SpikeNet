
function visualize_grid_firing_centre(R)
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

figure('Name','Vis','color','w','NumberTitle','off');

axis equal;
box on;
set(gca,'xtick',[],'ytick',[]);

xlim([-hw hw]);
ylim([-hw hw]);
hold on;

j = 1;

ang=0:0.01:2*pi;
x_shift_vs = [0 fw fw -fw -fw fw -fw 0 0 ];
y_shift_vs = [0 fw -fw fw -fw 0  0   fw -fw];

for t = 1:length(t_mid);
    h2 = plot(100,0);
    h3 = plot(100,0);
    
    ind_range_tmp = ind_a_vec(t):ind_b_vec(t);
    h1 = plot(x_pos_o(ind_range_tmp), y_pos_o(ind_range_tmp), 'bo');
    if sum(t_mid_chosen == t_mid(t)) == 1
        
        x_tmp = x_mean_chosen(j);
        y_tmp = y_mean_chosen(j);
        r_cos = x_tmp+dist_std_chosen(j)*cos(ang);
        r_sin = y_tmp+dist_std_chosen(j)*sin(ang);
        
        h3 = plot(r_cos - x_shift_vs(1),r_sin - y_shift_vs(1),'r', ...
            r_cos - x_shift_vs(2),r_sin - y_shift_vs(2),'r',...
            r_cos - x_shift_vs(3),r_sin - y_shift_vs(3),'r',...
            r_cos - x_shift_vs(4),r_sin - y_shift_vs(4),'r',...
            r_cos - x_shift_vs(5),r_sin - y_shift_vs(5),'r', ...
            r_cos - x_shift_vs(6),r_sin - y_shift_vs(6),'r', ...
            r_cos - x_shift_vs(7),r_sin - y_shift_vs(7),'r', ...
            r_cos - x_shift_vs(8),r_sin - y_shift_vs(8),'r', ...
            r_cos - x_shift_vs(9),r_sin - y_shift_vs(9),'r');
        
        if j == 1
            h2 = plot( x_tmp, y_tmp, 'r>', 'MarkerSize', 8);
        else
            h2 = plot( x_tmp, y_tmp, 'rx', 'MarkerSize', 16);
            xs = sprintf('jump_size = %4.2f, radius = %4.2f', jump_dist(j-1), dist_std_chosen(j));
            xlabel(xs);
        end
        
        
        j = j + 1;
    end
    pause(0.05);
    delete(h1);
    delete(h2);
    delete(h3);
    ts = sprintf('time = %8.1f ms', t_mid(t)*0.1);
    title(ts);
    
    
end
end
