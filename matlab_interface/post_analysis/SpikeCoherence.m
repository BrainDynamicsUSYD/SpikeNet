function c = SpikeCoherence(a, b, tau)

a = movsum(a,tau,'Endpoints','discard');
a = full(a(1:tau:end));

b = movsum(b,tau,'Endpoints','discard');
b = full(b(1:tau:end));

c = sum(a.*b) / sqrt( sum(a)*sum(b) );
end