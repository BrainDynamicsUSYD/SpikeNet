
% test multiple color for one line

x = 0:.05:2*pi;
y = sin(x);
z = zeros(size(x));
col = x;  % This is the color, vary with x in this case.

surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
