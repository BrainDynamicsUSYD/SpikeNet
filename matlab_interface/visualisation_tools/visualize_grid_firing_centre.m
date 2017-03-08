
function visualize_grid_firing_centre(R,mode)
close all;
clc;
% mode = 'bayesian';

hw = 31;
fw = 2*hw+1;


t_mid = R.grid.t_mid;
ind_a_vec = R.grid.ind_ab(1,:);
ind_b_vec = R.grid.ind_ab(2,:);

switch mode
    case 'bayesian'
        x_centre = R.grid.bayes.centre(1,:);
        y_centre = R.grid.bayes.centre(2,:);
        width = R.grid.bayes.radius;
    case 'quick'
        x_centre = R.grid.quick.centre(1,:);
        y_centre = R.grid.quick.centre(2,:);
        width = R.grid.quick.radius;
end


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

ang=0:0.01:2*pi;
x_shift_vs = [0 fw fw -fw -fw fw -fw 0 0 ];
y_shift_vs = [0 fw -fw fw -fw 0  0   fw -fw];

for t = 1:length(t_mid)
    h2 = plot(100,0);
    h3 = plot(100,0);
    
    ind_range_tmp = ind_a_vec(t):ind_b_vec(t);
    h1 = plot(x_pos_o(ind_range_tmp), y_pos_o(ind_range_tmp), 'bo');
    if ~isnan(x_centre(t))
        x_tmp = x_centre(t);
        y_tmp = y_centre(t);
        r_cos = x_tmp+width(t)*cos(ang);
        r_sin = y_tmp+width(t)*sin(ang);
        h3 = plot(r_cos - x_shift_vs(1),r_sin - y_shift_vs(1),'r', ...
            r_cos - x_shift_vs(2),r_sin - y_shift_vs(2),'r',...
            r_cos - x_shift_vs(3),r_sin - y_shift_vs(3),'r',...
            r_cos - x_shift_vs(4),r_sin - y_shift_vs(4),'r',...
            r_cos - x_shift_vs(5),r_sin - y_shift_vs(5),'r', ...
            r_cos - x_shift_vs(6),r_sin - y_shift_vs(6),'r', ...
            r_cos - x_shift_vs(7),r_sin - y_shift_vs(7),'r', ...
            r_cos - x_shift_vs(8),r_sin - y_shift_vs(8),'r', ...
            r_cos - x_shift_vs(9),r_sin - y_shift_vs(9),'r');
        
        h2 = plot( x_tmp, y_tmp, 'r>', 'MarkerSize', 8);
        
    end
    pause(0.05);
    delete(h1);
    delete(h2);
    delete(h3);
    ts = sprintf('time = %8.1f ms', t_mid(t)*0.1);
    title(ts);
end

end
