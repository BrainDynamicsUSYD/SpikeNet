function writeExplVarHDF5(FID, varargin)

% check input
var_num = length(varargin)/2;
if mod(var_num,1) ~= 0
    disp('wrong PopPara input format!\n')
    disp('var_num: ');disp(var_num);
else
    
    para_str = [];
    for i = 1:var_num
        para_str = [para_str, varargin{i*2-1}, ',' , num2str(varargin{i*2}), ',']; %#ok<AGROW>
    end
    
    hdf5write(FID,'/config/explanatory_variables', para_str,'WriteMode','append');
end

end

