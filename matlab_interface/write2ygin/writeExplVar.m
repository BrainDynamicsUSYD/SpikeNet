function writeExplVar(FID, varargin)


for i = 1:length(varargin)/2
    % Read var_name, var_value
    var_name = varargin{2*i-1};
    var_value = varargin{2*i};

    % var-name should be valid matlab-type variable name!
    fprintf(FID, '%s\n', '> explanatory variable');
    fprintf(FID, '%s,\n', var_name);
    if isnumeric(var_value)
        fprintf(FID, '%.9f,\n\n', var_value); % numeric
    else
        % No comma!!
        var_value = strrep(var_value, ',', ' '); % replace comma with whitespace
        fprintf(FID, '%s,\n\n', var_value); % string
    end
end

end

