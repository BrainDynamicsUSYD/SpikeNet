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

N = [3969, 1000] % replace this!
hw = 31;
fw = hw*2+1;
        
if ~isempty(files)
    for fi = 1:length(files)
        % OutData{id_out}.file = files{id_out};
        [file_dir, file_name, ~] = fileparts(files{fi});
        fprintf('Current file is: %s\n', files{fi});
        tic;
        fprintf('Reading ygin_syn file...');
        R_EE = read_ygin_syn(files{fi},1,1);
        W_EE = sparse(R_EE{1}.I, R_EE{1}.J, R_EE{1}.K, N(1), N(1));
        clear R_EE;
        R_IE = read_ygin_syn(files{fi},1,2);
        W_IE = sparse(R_IE{1}.I, R_IE{1}.J, R_IE{1}.K, N(1), N(2));
        clear R_IE;
        R_EI = read_ygin_syn(files{fi},2,1);
        W_EI = sparse(R_EI{1}.I, R_EI{1}.J, R_EI{1}.K, N(2), N(1));
        clear R_IE;
        fprintf('done.\n');
        
        

        pbc = 1; % periodic boundary condition
        EE_dist = GridDM([fw fw 1], [fw fw 1], pbc); % a full N1-by-N1 matrix of the distance between each pair of neurons
        r = 7;% radius of the local field
        EE_local = EE_dist <= r; % symmetric matrix
        W_EE_local_sum = vec2mat(diag(EE_local*W_EE*EE_local), fw);
        W_EIE_local_sum = vec2mat(diag(EE_local*W_IE*W_EI*EE_local), fw);

        save([file_name, '_EE_EIE.mat'], 'W_EE_local_sum', 'W_EIE_local_sum','r');
    end
    
end

% mat_EE = (mat_EE - mean(mat_EE(:)))/std(mat_EE(:));%standard normalized the matrix
% mat_EIE = (mat_EIE - mean(mat_EIE(:)))/std(mat_EIE(:));
% mat_EIE = -1*mat_EIE;

end