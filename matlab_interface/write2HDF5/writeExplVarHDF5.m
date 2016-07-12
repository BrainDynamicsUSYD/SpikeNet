function writeExplVarHDF5(FID, varargin)


for i = 1:length(varargin)/2
    % Read var_name, var_value
    var_name = varargin{2*i-1};
    var_value = varargin{2*i};

    hdf5write(FID,['/config/explanatory_variable/',var_name],var_value,'WriteMode','append');
end

end

