function [A] = generate_graph_given_in_out_degree(degree_in, degree_out, no_self)

if nargin == 2
    no_self = 1; % no self-connection
end

N = length(degree_in);
stub_in = degree2stub( degree_in );
stub_out = degree2stub( degree_out );

if no_self == 1
    self_num = length(stub_in);
    while self_num > 0
       swap = randperm(length(stub_in), 2);
       stub_in_new = stub_in;
       stub_in_new(swap) = stub_in_new(fliplr(swap));
       self_num_new = sum(stub_in_new == stub_out);
       if self_num_new <  self_num
           stub_in(swap) = stub_in(fliplr(swap));
           self_num = self_num_new;
       end
    end
end

A = full(sparse(stub_out,stub_in,1,N,N)); % generate adjacency matrix

end