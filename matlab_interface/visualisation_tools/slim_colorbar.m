function [ ax_c ] = slim_colorbar(varargin)
% SLIM_COLORBAR Display color bar (color scale)
% Just like matlab colorbar but it keeps the size of the associated axes
% unchanged and allows user to specify the width of colorbar.
% 
% USAGE:
%  slim_colorbar
%  slim_colorbar(width)
%  slim_colorbar(width, ax_acc)
%  [ ax_c ] = slim_colorbar(...)
% 
% ARGUMENTS (optional):
%        width -- width of colorbar (% of figure width, suggest 0.01)
%        ax_acc -- the associated axes (the default is the current one)
% 
% OUTPUT: 
%        ax_c -- axes handle of the colorbar
% 
% EXAMPLE:
% {
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
%     slim_colorbar(0.01, ax_acc); 
%     % or slim_colorbar(0.01)
%     % or simply slim_colorbar;
% }
%
% Yifan Gu, Feb 2017
% yigu8115@gmail.com

width = 0.01;
ax_acc = gca;

if nargin >= 1
    width = varargin{1};
end
if nargin == 2
    ax_acc = varargin{2};
end
    
ax_acc_size = get(ax_acc, 'position');
ax_c = colorbar('peer',ax_acc);
set(ax_acc, 'Position', ax_acc_size);pause(1)
ax_c_size = get(ax_c,'Position');
ax_c_size(3)= width;
set(ax_c,'Position',ax_c_size)

end 