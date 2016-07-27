function multicolor_line(x, y, color, linewidth)

if nargin == 3
    linewidth = 2;
end


x = x(:)';
y = y(:)';
color = color(:)';

z = zeros(size(x));
if length(color) == 1;
    color = color*ones(size(x));
end

surface([x;x],[y;y],[z;z],[color;color],...
    'facecol','no',...
    'edgecol','interp','linew', linewidth);

% % use the following commands to control the color axis when the plot 
% is not done in one go.

% caxis([0,c_max]);

end





