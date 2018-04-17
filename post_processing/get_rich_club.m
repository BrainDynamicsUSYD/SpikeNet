function R = get_rich_club(R)

disp('Getting rich club...')

n_sub = 100;
n_sample = 2000;

L = lattice_nD(2, 31);
I_sample = find(  L(:,1).^2+L(:,2).^2 <= n_sample/pi );
N_sample = length(I_sample);


s = strsplit(R.stamp,'_');
FID = [s{1},'_in.h5'];
[syn_cell] = read_syn_h5(FID, 1, 1, 1);

I = syn_cell{1}.I;
J = syn_cell{1}.J;
K = syn_cell{1}.K;
clear syn_cell;
W = full(sparse(I,J,K, 63^2, 63^2));
W = W(I_sample, I_sample);
clear I J K;

% topological rich club
in_degree = sum(W > 0,1);
out_degree = sum(W > 0,2);
Rw_t = rich_club_bd( double(W > 0) );
%
Rw_t_sub = [];
for i = 1:n_sub
    i
    [cij,flag] = makerandCIJdegreesfixed(in_degree(:),out_degree(:));
    Rw_sub_tmp = rich_club_bd( double(cij > 0) );
    Rw_t_sub = [Rw_t_sub; Rw_sub_tmp(:)']; %#ok<AGROW>
end


% weight rich club
Rw_w = rich_club_wd( W );
[I,J,K] = find(W);
clear W;
%
Rw_w_sub = [];
for i = 1:n_sub
    i
    W_sub = full(sparse(I, J, K(randperm(length(K))), 63^2, 63^2));
    Rw_sub_tmp = rich_club_wd( W_sub );
    Rw_w_sub = [Rw_w_sub; Rw_sub_tmp(:)']; %#ok<AGROW>
end

R.rich_club.Rw_t = Rw_t;
R.rich_club.Rw_t_sub = mean(Rw_t_sub);
R.rich_club.Rw_w = Rw_w;
R.rich_club.Rw_w_sub = mean(Rw_w_sub);

end





