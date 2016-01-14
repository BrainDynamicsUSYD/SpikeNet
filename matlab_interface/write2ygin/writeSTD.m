function writeSTD(FID, pre_pop_ind, post_pop_ind, STD_step)
% add STD to connection
%     FID: file id for writing data
% pop_ind: 

pre_pop_ind = pre_pop_ind - 1; % from matlab to c++ index
post_pop_ind = post_pop_ind - 1;

fprintf(FID, '%s\n', '> INIT008');
fprintf(FID, '%d, %d, %d,', pre_pop_ind, post_pop_ind, STD_step);
fprintf(FID,'\n\n');
end

