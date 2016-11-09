function writeSpikeFileInputPopHDF5(FID, fname,pop_ind)
% add file name of DVS input event/spike data
% the file is assumed to be an HDF5 file with contents 
% x, y, t, pol as lists of event data
%     FID: file id for writing data

pop_ind = pop_ind - 1; % from matlab to c++ index

% hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/file_spike_input/fname'],fname,'WriteMode','append');
% the hdf5write function seems incompatible with variable length strings
% which are used to read/write in the C++ simulator code, so the below is a
% much more complicated series of steps using low level functions to write
% the string to the hdf5 file as a variabel length string

% Create a new file using the default properties.
%
file = H5F.open (FID, 'H5F_ACC_RDWR', 'H5P_DEFAULT');%, 'H5P_DEFAULT');

%
% Create file and memory datatypes.  For this example we will save
% the strings as FORTRAN strings.
%
filetype = H5T.copy ('H5T_C_S1');
H5T.set_size (filetype,'H5T_VARIABLE');
memtype = H5T.copy ('H5T_C_S1');
H5T.set_size (memtype, 'H5T_VARIABLE');

%
% Create dataspace.  Setting maximum size to [] sets the maximum
% size to be the current size.
%
space = H5S.create_simple (1,fliplr( 1), []);

% Create group creation property list and set it to allow creation
% of intermediate groups.
%
gcpl = H5P.create ('H5P_LINK_CREATE');
H5P.set_create_intermediate_group (gcpl, 1);
% Create the group /G1/G2/G3.  Note that /G1 and /G1/G2 do not
% exist yet.  This call would cause an error if we did not use the
% previously created property list.
%
group = H5G.create (file, ['/config/pops/pop',num2str(pop_ind),'/file_spike_input'], gcpl,'H5P_DEFAULT', 'H5P_DEFAULT');


%
% Create the dataset and write the variable-length string data to
% it.
%
dset = H5D.create (file, ['/config/pops/pop',num2str(pop_ind),'/file_spike_input/fname'], filetype, space, 'H5P_DEFAULT');
H5D.write (dset, memtype, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', {fname});

%
% Close and release resources.
%
H5D.close (dset);
H5S.close (space);
H5T.close (filetype);
H5T.close (memtype);
H5F.close (file);


end

