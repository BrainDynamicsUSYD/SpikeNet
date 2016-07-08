function local_field_connections(varargin)
tic;
d_ygin_syn= dir('*ygin_syn');
hw = 31;% half width of the grid
fw = 2*hw + 1;
r = 7;%radius of the local field
mat_c_k = zeros(fw);
mat_EE = zeros(fw);
mat_EIE = zeros(fw);
[Lattice, N] = lattice_nD(2, hw);
loop_number = 7;
R_EE = read_ygin_syn(d_ygin_syn(loop_number).name,1,1);
R_IE = read_ygin_syn(d_ygin_syn(loop_number).name,1,2);
R_EI = read_ygin_syn(d_ygin_syn(loop_number).name,2,1);
R_EE = R_EE{1};
R_IE = R_IE{1};
R_EI = R_EI{1};
type1 = 'EE';
type2 = 'EIE';
for i = 1:N
    dist = lattice_nD_find_dist(Lattice, hw, i);
    local = (find(dist <= r));
    con_EE = ismember(R_EE.I,local).*ismember(R_EE.J,local);% EE
    con_EIE_IE = ismember(R_IE.I,local);%EIE_IE
    con_EIE_EI = ismember(R_EI.I,R_IE.J(con_EIE_IE));%EIE_EI
    con_EIE = ismember(R_EI.J(con_EIE_EI),local);%EIE
    mat_EE(i) = sum(R_EE.K(logical(con_EE)));
    mat_EIE(i) = mean(R_IE.K)*sum(R_EI.K(con_EIE));
end
%minmax EE[6.6874 7.1309] EIE[0.1426 0.1462] EIIE[0.0012 0.0012]
mat_EE = (mat_EE - mean(mat_EE(:)))/std(mat_EE(:));%standard normalized the matrix
mat_EIE = (mat_EIE - mean(mat_EIE(:)))/std(mat_EIE(:));
mat_EIE = -1*mat_EIE;
fprintf('Work before chosing the coefficients done...\n');
toc;
for k = 0.5
    mat_c_k = (1 - k)*mat_EE + k*mat_EIE;
    mat_c_k = permute(mat_c_k,[2 1]);
    mat_c_k = flipud(mat_c_k);
    mat_c_k = (mat_c_k - mean(mat_c_k(:)))/std(mat_c_k(:));%standard normalized the matrix
    figure(1)
    imagesc(mat_c_k)
    colorbar
    set(gcf,'renderer','zbuffer');
    t2 = sprintf(['local field connections %0.2g%s-%0.2g%s,loop number = %04i,r = %d'],1 - k,type1,k,type2,loop_number,r);
    title(t2);
    name_connections = sprintf(['%04i_local_field_connections_%0.2g%s-%0.2g%s_r%d.pdf'],loop_number,1 - k,type1,k,type2,r);
    saveas(gca,name_connections);
    fprintf('Saving results...\n');
end
fprintf('Done...\n');
toc;
end