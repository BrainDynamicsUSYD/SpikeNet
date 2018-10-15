function R = get_susceptibility( R ) 
%  Griffiths phase and long-range correlations in a biologically motivated
%  visual cortex model 

rho = R.reduced.num_spikes{1}/R.N(1);

susceptibility = R.N(1) * ( mean(rho.^2) - mean(rho)^2) / mean(rho);

% output resultss
R.susceptibility = susceptibility;

end


