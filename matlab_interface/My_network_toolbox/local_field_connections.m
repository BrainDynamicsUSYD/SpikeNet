function local_field_connections(files)
% Prepare filenames
if nargin == 0
    % If given no argument, search for matches under CURRENT directory
    dir_strut = dir('*.ygin_syn');
    num_files = length(dir_strut);
    files = cell(1,num_files);
    for fi = 1:num_files
        files{fi} = dir_strut(fi).name;
    end
    
end

if ~isempty(files)
    
    hw = 31;% half width of the grid
    fw = 2*hw + 1;
    [Lattice, N] = lattice_nD(2, hw);
    for fi = 1:length(files)
        % OutData{id_out}.file = files{id_out};
        [file_dir, file_name, ~] = fileparts(files{fi});
        
        fprintf('Current file is: %s\n', files{fi});
        
       
        r = 7;% radius of the local field
        
        sum_K_EE_local = zeros(fw);
        sum_K_EIE_local = zeros(fw);
        tic;
        fprintf('Reading ygin_syn file...');
        R_EE = read_ygin_syn(files{fi},1,1);
        R_IE = read_ygin_syn(files{fi},1,2);
        R_EI = read_ygin_syn(files{fi},2,1);
        fprintf('done.\n');
        toc;
        R_EE = R_EE{1};
        R_IE = R_IE{1};
        R_EI = R_EI{1};
        for ind = 1:N
            dist = lattice_nD_find_dist(Lattice, hw, ind);
            ind_E_local = (find(dist <= r)); % the indices of the neurons in pop E that are considered as local
            clear dist;
            
            i_EE_local = ismember(R_EE.I, ind_E_local); % true if the connection from pop E to pop E originates locally
            j_EE_local = ismember(R_EE.J, ind_E_local); % true if the connection from pop E to pop E terminates locally
            sum_K_EE_local(ind) = sum( R_EE.K(i_EE_local & j_EE_local) ); % the sum of the connection stengths from local neurons to local neurons in pop E
            clear i_EE_local j_EE_local;
             
            i_IE_local = ismember(R_IE.I, ind_E_local);% true if the connection from pop E to pop I originates locally in pop E
            j_EI_local = ismember(R_EI.J, ind_E_local);% true if the connection from pop I to pop E terminates locally in pop E
            
            ind_I_recv_local = unique(R_IE.J(i_IE_local)); % the indices of the neurons in pop I that receive connections from local neurons in pop E  
            ind_I_send_local = unique(R_EI.I(j_EI_local)); % the indices of the neurons in pop I that send connections to local neurons in pop E  
            
            ind_I_local = intersect(ind_I_recv_local, ind_I_send_local); % the indices of the neurons in pop I that are considered as local, which is defined as the intersection of the two above sets  
            clear ind_I_recv_local ind_I_send_local;
            
            sum_K_IE_local_EIE = sum( R_IE.K( i_IE_local & ismember(R_IE.J, ind_I_local) ) ); % the sum of the connection stengths from local neurons in pop E to local neurons in pop I
            sum_K_EI_local_EIE = sum( R_EI.K( j_EI_local & ismember(R_EI.I, ind_I_local) ) ); % the sum of the connection stengths from local neurons in pop I to local neurons in pop E
            clear ind_I_local i_IE_local j_EI_local;
            
            sum_K_EIE_local(ind) = sum_K_EI_local_EIE*sum_K_IE_local_EIE;
            
        end
        save([file_dir, '\', file_name, '_EE_EIE.mat'], 'sum_K_EE_local', 'sum_K_EIE_local','r');
    end
    
end

% mat_EE = (mat_EE - mean(mat_EE(:)))/std(mat_EE(:));%standard normalized the matrix
% mat_EIE = (mat_EIE - mean(mat_EIE(:)))/std(mat_EIE(:));
% mat_EIE = -1*mat_EIE;

end