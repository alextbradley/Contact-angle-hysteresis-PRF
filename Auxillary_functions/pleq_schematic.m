function pleq_schematic(dT, nu, xl)
%plot the equilibrium in schematic form

pleq(dT, nu, xl); %run the equilibrium plot

%remove any labels
xticks([]);
yticks([]);
box off
%ax = gca; ax.Visible = 'off'; %make invisible
ax = gca;
ax.XAxis.Visible = 'off'; %make X axis invisible

%add a solid line for the base
plot([0,0], [-1.1, 1.1], 'k', 'linewidth', 4)

%adjust ylim
ylim([-1.1, 1.1])

%make background white
set(gcf, 'color', 'w')

