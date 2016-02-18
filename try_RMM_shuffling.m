


A = transpose(R.spike_hist{1}(1:40,:));


% Now shuffle
% Implementation of the Raster Marginals Model, as introduced in
% "Population rate dynamics and multineuron firing patterns in 
%  sensory cortex", Journal of Neuroscience.

neuron_num = size(A,2); % keep it ~100
shuffle_num = 10*nchoosek(neuron_num,2); 
c = ceil(rand(shuffle_num,2)*neuron_num); % two randomly selected columns

for i = 1:shuffle_num

  I = A(:,c(i,1)) + A(:,c(i,2)) == 1; % where the 2 columns don't coincide
  cA = A(I,[c(i,1) c(i,2)]); % a copy of the part that matters, to make it run faster
  i01 = find(cA(:,1) == 0);
  i10 = find(cA(:,1) == 1);  
    
  toFlip = ceil(min(length(i01), length(i10))/2); % how many 01s & 10s to flip
  
  i01 = i01(randperm(length(i01)));
  i01 = i01(1:toFlip);
  i10 = i10(randperm(length(i10)));
  i10 = i10(1:toFlip);
  
  % the flip itself:
  cA(i01,1) = true; cA(i01,2) = false;
  cA(i10,1) = false; cA(i10,2) = true;
  A(I,[c(i,1) c(i,2)]) = cA;  

end;
A = A';



