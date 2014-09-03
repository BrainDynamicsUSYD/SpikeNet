function plot_comment(R)

comments = R.comments;
set(gca,'visible','off')
text(0.5, 0.5, comments, ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'center',...
    'FontSize',10,'FontWeight','normal', 'interpreter', 'none'); % ...'interpreter', 'none'... to show underscore

end