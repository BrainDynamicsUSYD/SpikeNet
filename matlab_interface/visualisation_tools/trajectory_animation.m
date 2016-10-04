function trajectory_animation(x,y, t_p, refresh_steps)
%TRAJECTORY_ANIMATION Generate animation of a 2-D time series data
%
% ARGUMENTS:
%        x -- A vector that contains the x data
%        y -- A vector that contains the y data
%        t_p -- the pause duration for each step (sec)
%        refresh_steps -- (Optional) the number of steps after which the
%           drawn trejectory is erased for performance concerns
%
% USAGE:
%{
% % Make up some data.
% x = 0.1*rand(1,3*10^3) + sin(linspace(0, 40*pi, 3*10^3));
% y = 0.1*rand(1,3*10^3) + cos(linspace(0, 40*pi, 3*10^3));
% refresh_steps = 700;
% t_p = 0.02;
% trajectory_animation(x,y, t_p, refresh_steps)
%}
%
% Yifan Gu, Sep 2016
% yigu8115@gmail.com

figure('NumberTitle','off','name','trajectory','color','w');
box on;
xlim(minmax(x))
ylim(minmax(y))
hold on;
ball = plot(x(1), y(1),'ro');
refresh_acc = 0;
if nargin < 3
    t_p = 0.01;
end
if nargin < 4
    refresh_steps = 500;
end
for i = 2:length(x)
    if refresh_acc == refresh_steps
        cla;
        ball = plot(x(i-1), y(i-1),'ro');
        refresh_acc = 0;
    end
    
    refresh_acc = refresh_acc + 1;
    plot([x(i-1) x(i)], [y(i-1) y(i)],'b');
    set(ball,'XData',x(i),'Ydata', y(i));
    xlabel(sprintf('step = %g',i))
    pause(t_p);
end

end