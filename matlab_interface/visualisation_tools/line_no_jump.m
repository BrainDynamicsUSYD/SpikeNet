function line_no_jump(x, y, marker, thre)
% thre is in percentage of the range of y
hold on;

if nargin < 4
    thre = 0.6;
end

y_diff = abs(diff(y));
jumps = [false y_diff > thre*range(y) true];

a = 1;
for i = 1:length(jumps)
    if jumps(i)
        b = i-1;
        plot(x(a:b), y(a:b), marker);
        a = i;
    end
end
    
end
