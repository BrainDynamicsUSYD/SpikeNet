function writeChemicalConnectionHDF5(FID, type, i_pre, j_post, I, J, K, D)
%      FID: file id for writing data
%     type: type of chemical connection (1:AMPA, 2:GABA, 3:NMDA)
%    i_pre: pre-synaptic population index
%   j_post: post-synaptic population index
%        I: pre-synaptic neuron index vector
%        J: post-synaptic neuron index vector
%        K: synaptic connection strength vector
%        D: synaptic delay vector
%
% Note that (I,J,K) defines a (sparse) adjacency matrix
%  and that (I,J,D) defines a delay matrix

% for C/C++ index convetion
type = type -1;
i_pre = i_pre -1;
j_post = j_post -1;
I = I-1; 
J = J-1;
if nargin == 7
    D = zeros(size(I)); % auto zero delay
end
% write (note: no white space!!!)
%fprintf(FID, '%s\n', '# chemical connection // (type,i_pre,j_post);I;J;K;D;');
%work out how many synapses already added
info=h5info(FID,'/config/syns');
n_syn=length(info.Groups);
n_syn=n_syn+1;
if n_syn==1
    h5create(FID,'/config/syns/n_syns',1,'Datatype','int32');
end
h5write(FID,'/config/syns/n_syns',int32(n_syn));

hdf5write(FID,['/config/syns/syn',num2str(n_syn),'/INIT006/type'],type,'WriteMode','append');
hdf5write(FID,['/config/syns/syn',num2str(n_syn),'/INIT006/i_pre'],i_pre,'WriteMode','append');
hdf5write(FID,['/config/syns/syn',num2str(n_syn),'/INIT006/j_post'],j_post,'WriteMode','append');
hdf5write(FID,['/config/syns/syn',num2str(n_syn),'/INIT006/I'],I,'WriteMode','append');
hdf5write(FID,['/config/syns/syn',num2str(n_syn),'/INIT006/J'],J,'WriteMode','append');
hdf5write(FID,['/config/syns/syn',num2str(n_syn),'/INIT006/K'],K,'WriteMode','append');
hdf5write(FID,['/config/syns/syn',num2str(n_syn),'/INIT006/D'],D,'WriteMode','append');

end