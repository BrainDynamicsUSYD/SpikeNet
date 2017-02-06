
function r = try_h5read(file, dataset_name)
r = [];
try
    r = h5read(file, dataset_name);
catch ME
    switch ME.identifier
        case 'MATLAB:imagesci:h5read:libraryError'
        otherwise
            rethrow(ME)
    end
    
end
end