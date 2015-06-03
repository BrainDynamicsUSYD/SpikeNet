function A=jacobi(f,x, delta)
% Computes the derivative (Jacobi matrix) A=Df(x) of the function f at x.
% 
%   Syntax: 
%           A = jacobi(f,x)
%   Arguments: 
%           f - function handle pointing to a function file  
%               which returns an mx1 matrix y=f(x).  
%           x - an nx1 matrix
%   Returns: 
%           A - an mxn matrix, the Jacobian of f 
%   Description:
%           The program jacobi uses a central difference quotient
%           to compute an approximation of the Jacobi matrix A=Df(x)
%           The file f.m must return an mx1 matrix y=f(x).  
%           The stepsize is h=1e-5.  
%           The optimal stepsize for the central difference quotient is 
%           h = (df/Kf)^(1/3), where df=1e-16 is the error in f(x) and 
%           Kf is the constant in Taylor's formula,
%           f(x)=f(a)+f'(a)(x-a)+0.5 f''(a)(x-a)^2 + Ef, |Ef| < Kf|x-a|^3.
%           The error is then 2(df)^(2/3)(Kf)^(1/3), approximately 1e-10.
%   Examples: 
%          z = jacobi(@sin, 0) returns 1
%          z = jacobi(@funk1, 3) returns 6
%          A = jacobi(@funk2, [1;2;3]) returns [2 4 6]


%----------------------------------------------------------------------
%   Date created: 1999-12-05,  updated: 2006-11-07
%   
%   Author:  Stig Larsson
%------------------------------+---------------------------------------+

% delta = 1e-5;   % stepsize in the central difference approximation

% test if f is a column vector
y=f(x);
if size(y,2)~=1 
  disp('jacobi: f is not a column vector'); 
  return;
end

n=length(x);

for j = 1:n

   e = zeros(n,1);
   e(j) = 1;                 % e is the j:th basis vector
   h = delta*e;              % h is a step in the j:th direction

   % compute the j:th column of A 
   A(:,j) = ( f(x+h) - f(x-h) )/(2*delta);

end

