function multicolor_line(x,y,color, linewidth)

if nargin == 3
    linewidth = 2;
end
x = x(:)';
y = y(:)';

z = zeros(size(x));
if length(color) == 1;
    color = color*ones(size(x));
end

surface([x;x],[y;y],[z;z],[color;color],...
    'facecol','no',...
    'edgecol','interp','linew', linewidth);
caxis([0,1]);
colormap('jet');

end





