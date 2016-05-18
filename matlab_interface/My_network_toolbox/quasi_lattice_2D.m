function Lattice = quasi_lattice_2D(N, hw)
% N is the number of points
% hw is the half-width, or 2*hw+1 is the full width
% The quasi-lattice is centred at (0,0)
% It is quasi because 
% (1) the coordinates of the points are uniformly distributed random numbers
% (2) the indices of the points are re-arranged to resemble a regular square
% lattice for convenience.

w = 2*hw + 1;

X = zeros(N,1);
while X(1) == 0 || X(end) == w
    X = sort(rand(N,1) * w);
end

segs = linspace(0, w, round(sqrt(N))+1);
Y = [];
for i = 1:round(sqrt(N))
    X_seg = X(X >= segs(i) & X < segs(i+1));
    Y_tmp = rand(length(X_seg),1) * w;
    [Y_tmp, Y_tmp_ind] = sort(Y_tmp);
    X(X >= segs(i) & X < segs(i+1)) = X_seg(Y_tmp_ind);
    Y = [Y; Y_tmp]; %#ok<AGROW>
end

Lattice = [X Y] - (w/2); % centred at (0,0)


% % visualization
% figure(1); axis([-hw hw -hw hw]);
% hold on
% for i = 1:N
%     plot(Lattice(i,1), Lattice(i,2),'o');pause(0.01);
% end


end

