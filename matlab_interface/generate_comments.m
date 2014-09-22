
function [R] = generate_comments(R)

fprintf('\t Gennerating auto comments...\n');
% comment line length
cmt_llength = 100;
% dump field
ExplVar = R.ExplVar;
step_killed = R.step_killed;
Hz_overall = R.Analysis.Hz_overall;
% Comments
comments = ' ';
if step_killed >= 0
comments = sprintf('Runaway killed at step %g, ', step_killed);
end
comments = [comments sprintf('Mean firing rate %.4g Hz, ', Hz_overall)];
% Expl variables
fname_cell = fieldnames(ExplVar);
for f = 1:length(fname_cell)
    fn = fname_cell{f};
    if isnumeric(ExplVar.(fn))
        comments = [comments sprintf('%s = %.4g, ', fn, ExplVar.(fn)) ];
    else
        comments = [comments sprintf('%s = %s, ', fn, ExplVar.(fn)) ];
    end
end
% break comments into multiple lines if too long
line_num = ceil(length(comments) / cmt_llength);
if line_num > 1
    comments = [comments repmat(' ', 1, line_num*cmt_llength-length(comments)) ]; % pad with white space
    comment_lines = cell(0,1);
    for ll = 1:line_num
        comment_lines{ll} = comments((ll-1)*cmt_llength+1:ll*cmt_llength);
    end
    comments = comment_lines;
end
% record results
R.comments = comments;

end

