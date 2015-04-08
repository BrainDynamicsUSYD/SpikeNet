


function [out]=get_eq3_4(miu_V, sigma_V, arg)
% numerical solution to eq(3.4)
% in the paper Mean-driven and fluctuation-driven persistent activity in
% recurrent networks

t_m = 10; %ms
t_ref = 2; %ms

V_th = 20; %mV
V_res = 10; %mV is V_res the same as V_r ?? Yes.


upper = (V_th-miu_V)/(sqrt(2)*sigma_V);
lower = (V_res-miu_V)/(sqrt(2)*sigma_V);


% % symbolic functions are not very reliable !!!
% syms x y;
% 
% f1 = exp(x^2)*(1+erf(x));
% v = 1/(t_ref + t_m*sqrt(pi)*double(int(f1, x, lower, upper)));  % double() does numerical integration
% % the above result is correct and the unit is KHz = (ms)^-1
% % error message says try vpa instead of double(); what's vpa()??
% f2 = exp(y^2)*((1+erf(y))^2);
% f3 = exp(x^2)*int(f2, y, -inf, x);
% CV2 = 2*pi*v^2*double( int(f3, x, lower, upper)); % the value of it is not right!!!!! why?!



%% another approach: 2D integral over some non-rectangle area, also faster
fx = @(x) exp(x.^2).*(1+erf(x));
v = 1/(t_ref + t_m*sqrt(pi)*integral(fx, lower, upper)); 
fxy = @(x,y) (y<x).*exp(x.^2).*exp(y.^2).*((1+erf(y)).^2);
CV2 = 2*pi*v^2*integral2(fxy , lower, upper, -10,  upper); % -3 instead of -inf for numerical approximation
% be careful about the arbitrariness of the value -3!!!!!


%%
CV2 = CV2*t_m^2; % my guess, which took me 3hr of reading and worked!

CV = CV2^0.5;

if strcmp(arg,'v')
    out = v;
elseif strcmp(arg,'CV')
    out = CV;
elseif strcmp(arg,'both')
    out = [v, CV];
end

end
