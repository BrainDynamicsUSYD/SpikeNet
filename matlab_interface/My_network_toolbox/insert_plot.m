function [insert_axes_new] = insert_plot(main_fig, main_axes, insert_axes, x_insert, y_insert, insert_size)
%
% main_axes could be a subplot of main_fig
% otherwise use main_axes = findobj(main_fig,'Type','axes') or main_axes =
% get(gca)

% Defualt position parameters
if nargin == 3
    x_insert = 0.6; % [0,1]
    y_insert = 0.1; % [0,1]
end
if nargin < 6
    insert_size=0.35; % [0,1]
end


% Copy axes
insert_axes_new = copyobj(insert_axes,main_fig); 
% new_axes_handle = copyobj(old_axes_handle, fig_handle).
% the only difference between new_axes_handle and old_axes_handle is that the
% parent changes into fig_handle
% Note that axes object is child of figure object.

% set axes position
position_main = get(main_axes,'Position');
position_insert = get(insert_axes_new, 'Position');

left = position_main(1)+x_insert*position_main(3);
bottom = position_main(2)+y_insert*position_main(4);
width = position_main(3)*insert_size;
height = position_main(4)*insert_size;

set(insert_axes_new,'Position', [left, bottom, width, height]) %, 'FontUnits', 'normalized', 'FontSize',1/10);

delete(insert_axes);
end


