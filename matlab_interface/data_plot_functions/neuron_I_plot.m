function neuron_I_plot(R, pop_ind, sample_ind, current_type)
% current_type can be any non-empty subset of 
% {'I_leak','I_AMPA','I_GABA','I_NMDA','I_GJ','I_ext'}

marker = {'r','g','b','y','k','c','m'}; % for 6 different currents

% Dump fields
dt = R.dt;
step_tot = R.step_tot;
T = (1:step_tot)*dt;

% plot
hold on;
for type = 1:length(current_type)
        I = R.neuron_sample.(current_type{type}){pop_ind}(sample_ind,:);
        plot(T,I,marker{type});
end
lh = legend(current_type);
set(lh,'Interpreter','none'); % plain text instead of Tex

set(gca,'box','off')

end