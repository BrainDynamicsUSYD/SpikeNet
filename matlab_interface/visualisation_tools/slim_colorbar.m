function [ ax_c ] = slim_colorbar(ax_acc, width)
%SLIM_COLORBAR Display color bar (color scale)
% Just like matlab colorbar but it keeps the size of the associated axes
% unchanged and allows user to specify the width of colorbar.
%
% ARGUMENTS:
%        ax_acc -- the associated axes
%        width -- width of colorbar (% of figure width, suggest 0.01)
%
% OUTPUT: 
%        ax_c -- axes handle of the colorbar
% USAGE:
%{
%     figure('NumberTitle','Off','Name','Slim colorbar','color','w');
%     subplot(2,2,1);
%     plot(1:10);
%     subplot(2,2,3);
%     plot(1:10);
%     xlabel('Matlab colorbar')
%     colorbar
%     subplot(2,2,2);
%     plot(1:10);
%     ax_acc = subplot(2,2,4);
%     plot(1:10);
%     xlabel('Slim colorbar')
%     slim_colorbar(ax_acc, 0.01)
%}
%
% Yifan Gu, Sep 2016
% yigu8115@gmail.com

ax_acc_size = get(ax_acc, 'position');
ax_c = colorbar('peer',ax_acc);
set(ax_acc, 'Position', ax_acc_size);pause(1)
ax_c_size = get(ax_c,'Position');
ax_c_size(3)= width;
set(ax_c,'Position',ax_c_size)

end 