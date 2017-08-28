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
modify = 0;
writeBasicParaHDF5(FID, dt, step_tot, N, modify);


% initial condition settings
exteral_init_V = randn(1,N(1))*5 - 65;
writeExtInitVHDF5(FID, 1, exteral_init_V);

exteral_init_V = randn(1,N(2))*5 - 60;
writeExtInitVHDF5(FID, 2, exteral_init_V)

end
