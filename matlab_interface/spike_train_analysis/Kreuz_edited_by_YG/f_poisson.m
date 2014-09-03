function  poiss = f_poisson(len,rate,refrac)

uniform=rand(1,len);
poiss=refrac-log(1-uniform)/rate;