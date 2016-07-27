function add_figure_label( label, anno_shift, anno_fontsize )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% annotation
p = get(gca,'position');
p_anno = [p(1),p(2)+p(4), 0.01,0.01] + [anno_shift(:)', 0, 0];
annotation('textbox',p_anno,'string',label,'fontsize',anno_fontsize,'fontweight','bold','edgecolor','none');


end

