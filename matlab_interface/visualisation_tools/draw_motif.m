function draw_motif( centre, r, ind, node_size, arrow_size)
% see Motifs in Brain Networks, PLoS Biology 2004

if nargin < 4
node_size = 15;
arrow_size = 7.5;
end

hold on;

%       o 1
%
%
%
% o 2          o 3
%

lim = axis;

pa = get(gca,'position');
pf = get(gcf,'position');
y_scale_axis = (pa(4)*pf(4)) / (pf(3)*pa(3));
y_scale_coord = diff(ylim)/diff(xlim);
y_scale = y_scale_coord/y_scale_axis;


link_mat = [1 2;
    1 3;
    2 1;
    2 3;
    3 1;
    3 2];
ind2link_mat = {[3 5], [3 6], [3 4], [1 3 5], [3 5 6],...
    [1 2 3], [2 3 6], [1 3 5 6], [1 2 3 5],...
    [1 3 4 5], [3 4 5 6], [1 2 3 5 6], [1 2 3 4 5 6]};

for j = ind
    
    x = centre(j, 1) + r*cosd([90 210 -30]);
    y = centre(j, 2) + r*y_scale*sind([90 210 -30]);
    plot(x,y,'k.','MarkerSize', node_size);
    
    links = link_mat(ind2link_mat{j}, :);
    for i = 1:length(links(:,1))
        arrow([x(links(i,1)) y(links(i,1))],[x(links(i,2)) y(links(i,2))],arrow_size,'Color','k');
    end
end

axis(lim)

end