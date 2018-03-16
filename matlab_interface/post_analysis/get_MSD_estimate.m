function MSD = get_MSD_estimate(x,y,steps)


L = length(x);
N = floor(L/steps);
MSD = [];
for i = 1:N
    a = (i-1)*steps+1;
    b = (i-1)*steps+steps;
    x_tmp = x(a:b);
    y_tmp = y(a:b);
    D = (x_tmp-x_tmp(1)).^2 + (y_tmp-y_tmp(1)).^2;
    MSD = [MSD; D(:)'];
end

MSD = nanmean(MSD);
end