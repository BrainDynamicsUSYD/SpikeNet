function [Hal,He,Ht,pval95] = hurst(x,d,fontsize)
%HURST Calculate the Hurst exponent using R/S analysis.
%   H = HURST(X) calculates the Hurst exponent of time series X using 
%   the R/S analysis of Hurst [2], corrected for small sample bias [1,3,4]. 
%   If a vector of increasing natural numbers is given as the second input 
%   parameter, i.e. HURST(X,D), then it defines the box sizes that the 
%   sample is divided into (the values in D have to be divisors of the 
%   length of series X). If D is a scalar (default value D = 50) it is 
%   treated as the smallest box size that the sample can be divided into. 
%   In this case the optimal sample size OptN and the vector of divisors 
%   for this size are automatically computed. 
%   OptN is defined as the length that possesses the most divisors among 
%   series shorter than X by no more than 1%. The input series X is 
%   truncated at the OptN-th value. 
%   [H,HE,HT] = HURST(X) returns the uncorrected empirical and theoretical 
%   Hurst exponents.
%   [H,HE,HT,PV95] = HURST(X) returns the empirical 95% confidence 
%   intervals PV95 (see [4]).
%
%   If there are no output parameters, the R/S statistics is automatically 
%   plotted against the divisors on a loglog paper and the results of the 
%   analysis are displayed in the command window. HURST(X,D,FONTSIZE) 
%   allows to specify a fontsize different than 14 in the plotted figure.
%
%   References:
%   [1] A.A.Anis, E.H.Lloyd (1976) The expected value of the adjusted 
%   rescaled Hurst range of independent normal summands, Biometrica 63, 
%   283-298.
%   [2] H.E.Hurst (1951) Long-term storage capacity of reservoirs, 
%   Transactions of the American Society of Civil Engineers 116, 770-808.
%   [3] E.E.Peters (1994) Fractal Market Analysis, Wiley.
%   [4] R.Weron (2002) Estimating long range dependence: finite sample 
%   properties and confidence intervals, Physica A 312, 285-299.

%   Written by Rafal Weron (2011.09.30). 
%   Based on functions hurstal.m, hurstcal.m, finddiv.m, findndiv.m 
%   originally written by Witold Wiland & Rafal Weron (1997.06.30, 
%   2001.02.01, 2002.07.27).  

if nargin<3, 
    fontsize = 14; 
end
if nargin<2, 
    d = 50; 
end
if max(size(d)) == 1, 
    % For scalar d set dmin=d and find the 'optimal' vector d
    dmin = d;
    % Find such a natural number OptN that possesses the largest number of 
    % divisors among all natural numbers in the interval [0.99*N,N] 
    N = length(x); 
    N0 = floor(0.99*N);
    dv = zeros(N-N0+1,1);
    for i = N0:N,
        dv(i-N0+1) = length(divisors(i,dmin));
    end
    OptN = N0 + find(max(dv)==dv) - 1;
    % Use the first OptN values of x for further analysis
    x = x(1:OptN);
    % Find the divisors of x
    d = divisors(OptN,dmin);
else
    OptN = length(x);
end

N = length(d);
RSe = zeros(N,1);
ERS = zeros(N,1);

% Calculate empirical R/S
for i=1:N;
   RSe(i) = RScalc(x,d(i));
end

% Compute Anis-Lloyd [1] and Peters [3] corrected theoretical E(R/S)
% (see [4] for details)
for i=1:N;
    n = d(i); 
    K = [1:n-1];
    ratio = (n-0.5)/n * sum(sqrt((ones(1,n-1)*n-K)./K));
    if (n>340)
        ERS(i) = ratio/sqrt(0.5*pi*n);
    else
        ERS(i) = (gamma(0.5*(n-1))*ratio) / (gamma(0.5*n)*sqrt(pi));
    end
end

% Calculate the Anis-Lloyd/Peters corrected Hurst exponent
% Compute the Hurst exponent as the slope on a loglog scale
ERSal = sqrt(0.5*pi.*d);
Pal = polyfit(log10(d),log10( RSe - ERS + ERSal ),1);
Hal = Pal(1);
% Calculate the empirical and theoretical Hurst exponents
Pe = polyfit(log10(d),log10(RSe),1);
He = Pe(1);
P = polyfit(log10(d),log10(ERS),1);
Ht = P(1);

% Compute empirical confidence intervals (see [4])
L = log2(OptN);
% R/S-AL (min(divisor)>50) two-sided empirical confidence intervals
pval95 = [0.5-exp(-7.33*log(log(L))+4.21) exp(-7.20*log(log(L))+4.04)+0.5];
C = [   0.5-exp(-7.35*log(log(L))+4.06) exp(-7.07*log(log(L))+3.75)+0.5 .90];
C = [C; pval95                                                          .95];
C = [C; 0.5-exp(-7.19*log(log(L))+4.34) exp(-7.51*log(log(L))+4.58)+0.5 .99];

% Display and plot results if no output arguments are specified
if nargout < 1,
    % Display results
    disp('---------------------------------------------------------------')
    disp(['R/S-AL using ' num2str(length(d)) ' divisors (' num2str(d(1)) ',...,' num2str(d(length(d))) ...
        ') for a sample of ' num2str(OptN) ' values'])
    disp(['Corrected theoretical Hurst exponent    ' num2str(0.5,4)]);
    disp(['Corrected empirical Hurst exponent      ' num2str(Hal,4)]);
    disp(['Theoretical Hurst exponent              ' num2str(Ht,4)]);
    disp(['Empirical Hurst exponent                ' num2str(He,4)]);
    disp('---------------------------------------------------------------')

    % Display empirical confidence intervals
    disp('R/S-AL (min(divisor)>50) two-sided empirical confidence intervals')
    disp('--- conf_lo   conf_hi   level ---------------------------------')
    disp(C)
    disp('---------------------------------------------------------------')

    % Plot R/S
    h2 = plot(log10(d),log10(ERSal/(ERS(1)/RSe(1))),'b-');
    if fontsize > 10, 
        set(h2,'linewidth',2); 
    end;
    hold on
    h1 = plot(log10(d),log10(RSe-ERS+ERSal),'ro-');
    if fontsize > 10, 
        set(h1,'linewidth',2); 
    end;
    hold off
    set(gca,'Box','on','fontsize',fontsize);
    xlabel('log_{10}n','fontsize',fontsize);
    ylabel('log_{10}R/S','fontsize',fontsize);
    legend('Theoretical (R/S)','Empirical (R/S)')
end

function d = divisors(n,n0)
% Find all divisors of the natural number N greater or equal to N0
i = n0:floor(n/2);
d = find((n./i)==floor(n./i))' + n0 - 1;

function rs = RScalc(Z,n)
% Calculate (R/S)_n for given n
m = length(Z)/n;
Y = reshape(Z,n,m);
E = mean(Y);
S = std(Y);
for i=1:m;
    Y(:,i) = Y(:,i) - E(i);
end;
Y = cumsum(Y);
% Find the ranges of cummulative series
MM = max(Y) - min(Y);
% Rescale the ranges by the standard deviations
CS = MM./S;
rs = mean(CS);