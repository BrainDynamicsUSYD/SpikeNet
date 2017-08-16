function [e]=mlf(alf,bet,c,fi)
%
% MLF -- Mittag-Leffler function. 
%           MLF (alpha,beta,Z,P) is the Mittag-Leffler function E_{alpha,beta}(Z) 
%           evaluated with accuracy 10^(-P) for each element of Z. 
%           alpha and beta are scalars, P is integer, Z can be a vector or
%           a two-dimensional array. The output is of the same size as Z.
% (C) 2001-2012 Igor Podlubny, Martin Kacenak
% Last update: 2012-09-07


[cRows, cCols] = size(c);
c=m2v(c);

if nargin<4 , fi=6;		end
if nargin<3 || alf<=0 || fi<=0
else 
   [r,s]=size(c); [r1,s1]=size(alf); [r2,s2]=size(bet);
   mx=max([r,s]); mx1=max([r1,s1]); mx2=max([r2,s2]);
   if (r>1 && s>1) || (r1>1 && s1>1) || (r2>1 && s2>1) || (mx1>1 && mx2>1)
      sprintf('wrong number of input parameters')
   else
      if mx1>mx2 , mxx=mx1; e=zeros(mx,mx1);
      else mxx=mx2; e=zeros(mx,mx2);end;
      for i1= 1:mx
         for i2=1:mxx
            
            if r>s , z=c(i1,1); else z=c(1,i1); end
            if mx1>mx2 , if r1>s1 , alfa=alf(i2,1); else alfa=alf(1,i2);end, beta=bet;
            else if r2>s2 ,beta=bet(i2,1); else beta=bet(1,i2); end, alfa=alf; end
            if beta<0 , rc=(-2*log(10^(-fi)*pi/(6*(abs(beta)+2)*(2*abs(beta))^(abs(beta)))))^alfa;
               else  rc=(-2*log(10^(-fi)*pi/6))^alfa; end
            r0=max([1,2*abs(z),rc]);
            if (alfa==1 && beta==1)
               e(i1,i2)=exp(z);
            else
               if (alfa<1 && abs(z)<=1) || ( (1<=alfa && alfa <2) && abs(z)<=floor(20/(2.1-alfa)^(5.5-2*alfa))) || (alfa>=2 && abs(z)<=50)
                  oldsum=0;
                  k=0;
                  while (alfa*k+beta)<=0 
                     k=k+1;
                  end
                  newsum=z^k/gamma(alfa*k+beta);
                  while newsum~=oldsum
                     oldsum=newsum;
                     k=k+1;
                     term=z^k/gamma(alfa*k+beta);
                     newsum=newsum+term;
                     k=k+1;
                     term=z^k/gamma(alfa*k+beta);
                     newsum=newsum+term;
                  end
                  e(i1,i2)=newsum;
               else
                  if (alfa<=1 && abs(z)<=fix(5*alfa+10))
                     if ((abs(angle(z))>pi*alfa) && (abs(abs(angle(z))-(pi*alfa))>10^(-fi)))
                        if beta<=1
                           e(i1,i2)=rombint('K',0,r0,fi,alfa,beta,z);
                        else
                          eps=1;
                          e(i1,i2)=rombint('K',eps,r0,fi,alfa,beta,z)+ ...
                             rombint('P',-pi*alfa,pi*alfa,fi,alfa,beta,z,eps);
                       end
                    elseif (abs(angle(z))<pi*alfa && abs(abs(angle(z))-(pi*alfa))>10^(-fi))
                       if beta<=1
                           e(i1,i2)=rombint('K',0,r0,fi,alfa,beta,z)+ ...
                              (z^((1-beta)/alfa))*(exp(z^(1/alfa))/alfa);
                       else
                           eps=abs(z)/2;
                           e(i1,i2)=rombint('K',eps,r0,fi,alfa,beta,z)+ ...
                              rombint('P',-pi*alfa,pi*alfa,fi,alfa,beta,z,eps)+ ...
                              (z^((1-beta)/alfa))*(exp(z^(1/alfa))/alfa);
                        end
                     else
                        eps=abs(z)+0.5;
                        e(i1,i2)=rombint('K',eps,r0,fi,alfa,beta,z)+ ...
                           rombint('P',-pi*alfa,pi*alfa,fi,alfa,beta,z,eps);
                     end
                  else
                     if alfa<=1
                        if (abs(angle(z))<(pi*alfa/2+min(pi,pi*alfa))/2)
                           % alfa
                           newsum=(z^((1-beta)/alfa))*exp(z^(1/alfa))/alfa;
                           for k=1:floor(fi/log10(abs(z)))
                              newsum=newsum-((z^(-k))/gamma(beta-alfa*k));
                              % k
                           end
                           e(i1,i2)=newsum;
                        else
                           newsum=0;
                           for k=1:floor(fi/log10(abs(z)))
                              newsum=newsum-((z^-k)/gamma(beta-alfa*k));
                           end
                           e(i1,i2)=newsum;
                        end
                     else
                        if alfa>=2
                           m=floor(alfa/2);
                           sum=0;
                           for h=0:m
                              zn=(z^(1/(m+1)))*exp((2*pi*1i*h)/(m+1));
                              sum=sum+mlf(alfa/(m+1),beta,zn,fi);
                           end
                           e(i1,i2)=(1/(m+1))*sum;
                        else
                           e(i1,i2)=(mlf(alfa/2,beta,z^(1/2),fi)+mlf(alfa/2,beta,-z^(1/2),fi))/2;
                        end
                     end
                  end
               end
            end
         end
      end
   end
end

if isreal(c) 
    e = real(e);  
end

e = v2m(e,cRows,cCols);


function [res]=rombint(funfcn,a,b,order,varargin)
if nargin<4 ,order=6; end
if nargin<3
   Warning ('Error in input format')
else 
   rom=zeros(2,order);
   h=b-a;
   rom(1,1)=h*(feval(funfcn,a,varargin{:})+feval(funfcn,b,varargin{:}))/2;
   
   ipower=1;
   for i= 2:order
      sum=0;
      for j=1:ipower
      	sum=sum+feval(funfcn,(a+h*(j-0.5)),varargin{:});
      end
      rom(2,1)=(rom(1,1)+h*sum)/2;
      for k=1:i-1
         rom(2,k+1)=((4^k)*rom(2,k)-rom(1,k))/((4^k)-1);
      end
      
      for j=0:i-1
         rom(1,j+1)=rom(2,j+1);
      end
      ipower=ipower*2;
      h=h/2;
   end
   res=rom(1,order);
end


function res=K(r,alfa,beta,z)
res=r.^((1-beta)/alfa).*exp(-r.^(1/alfa)).*(r*sin(pi*(1-beta))-...
z*sin(pi*(1-beta+alfa)))/(pi*alfa*(r.^2-2*r*z*cos(pi*alfa)+z.^2));


function res=P(r,alfa,beta,z,eps)
w=(eps^(1/alfa))*sin(r/alfa)+r*(1+(1-beta)/alfa);
res=((eps^(1+(1-beta)/alfa))/(2*pi*alfa))*((exp((eps^(1/alfa))*cos(r/alfa)).*...
   (cos(w)+1i*sin(w))))/(eps*exp(1i*r)-z);

function A = v2m(V, M, N)
if numel(V)==M*N, 
    A = reshape(V, [N, M]);
    A = A' ;
else
    warning('Wrong dimensions of the output in V2M.')
end

function V = m2v(A)
M = A'; V = M(:);


