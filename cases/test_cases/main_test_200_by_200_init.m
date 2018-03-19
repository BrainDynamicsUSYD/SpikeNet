function main_test_200_by_200_init(varargin)

PBS_ARRAYID = varargin{1};
FID = dir([ sprintf('%04g-', PBS_ARRAYID),'*_in.h5']);
FID = FID.name;

dt = 0.1;
sec = round(10^3/dt); % 1*(10^3/dt) = 1 sec
step_tot = 1*sec; % use 10 second!

hw = 101;
% sptially embedded network
hw_0 = 31; % half-width, (31*2+1)^2 = 3969 ~ 4000, hw=44 gives 7921
N_e_0 = (hw_0*2+1)^2; %
N_i_0 = 1000;
%%%%%%%%
N_e = (hw*2+1)^2; %
N_i = round(N_i_0/N_e_0*N_e);
N = [N_e, N_i];


% write basic parameters
modify = 1;
writeBasicParaHDF5(FID, dt, step_tot, N, modify);


% initial condition settings
exteral_init_V = randn(1,N(1))*7 - 70;
exteral_init_V = ones(1,N(1))*(-70);
exteral_init_V(randperm(N(1), round(N(1)*0.05))) = 0;
writeExtInitVHDF5(FID, 1, exteral_init_V, modify);

exteral_init_V = randn(1,N(2))*7 - 70;
exteral_init_V = ones(1,N(2))*(-70);
exteral_init_V(randperm(N(2), round(N(2)*0.05))) = 0;
writeExtInitVHDF5(FID, 2, exteral_init_V, modify)

end
