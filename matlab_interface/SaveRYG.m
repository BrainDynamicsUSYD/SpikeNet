function SaveRYG( Result_cell )
%This function saves Result_cell as explanatory and response
%variables to .mat file
%   Detailed explanation goes here

Result_num = length(Result_cell);
for r_num = 1:Result_num
    R_temp = Result_cell{r_num};
    name_temp = strcat( Result_cell{r_num}.stamp, '_RYG.mat');
    fprintf('Saving results into file %s...\n', name_temp);
    save(name_temp, '-struct', 'R_temp', '-v7.3'); % -v7.3 for >2GB
end
disp('Saving done.')

end

