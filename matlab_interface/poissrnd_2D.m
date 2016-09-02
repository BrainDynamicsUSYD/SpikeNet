
function pois = poissrnd_2D(lambda, N, r)
%Yahav, Inbal, and Galit Shmueli. 
% "On generating multivariate Poisson data 
% in management science applications." 
% Applied Stochastic Models in Business and Industry 28.1 (2012): 91-102. APA	
%
% Yifan Gu, School of Physics, USYD, Aug 2016.
% yigu8115@gmail.com

D = 2;
if min(lambda) <= 5
    r = CorrectInitialCorrel(lambda(1), lambda(2), r);
end
var = [1 r; r 1];

max_count = 100;
count = 0;
err = 1;
err_max = 0.1;
while  err > err_max
    count = count + 1;
    if count > max_count
        break;
    end
    
    normal = mvnrnd(zeros(1, D), var, N);
    pois = normal;
    p = normcdf(normal);
    for s = 1:D
        pois(:,s) = poissinv( p(:,s), lambda(s) );
    end
    
    r_a = corrcoef(pois(:,1), pois(:,2));
    r_a = r_a(1,2);
    if r > 0.1
        err = abs(r_a-r)/r;
    else
        err = abs(r_a-r);
    end
    
end

if count > max_count
    pois = NaN;
    warning('Error cannot converge!')
end

end


% only an approximation
function corrected = CorrectInitialCorrel(lambda1, lambda2, r)
samples = 500;
u = rand(samples) ;

maxcor = corrcoef(poissinv(u, lambda1), poissinv(u, lambda2));
maxcor = maxcor(1,2);
mincor = corrcoef(poissinv(u, lambda1), poissinv(1-u, lambda2));
mincor = mincor(1,2);

a = -maxcor*mincor/(maxcor+mincor);
b = log( (maxcor+a)/a);

corrected = log((r+a)/a)/b;

if corrected > 1 || corrected < -1
    corrected = NaN;
end

end

