%Produce equilibria maps in (X+, \nu) space with colour according to
%\lambda_e for V = 0.2, 0.3, 0.4, 0.5 (each different volume on a
%different plot). Colours saturated at \lambda_e = 0.1 (note that,
%depending on matlab distribution, plots may have yellow patches in bottom
%left and right corners. these are artefacts of data collection technique:
%we do not go all the way down to zero.
%Note that although the plot only extends to those equilibria with
%\lambdae=0.3, analysis of similar plots with data up to large \lambda
%suggests that the no-equilibrium contour is approximately in the same
%place as at the border in this plot

figpref(4);
addpath('./Data')
fnames = ["V_pt2.mat","V_pt3.mat", "V_pt4.mat", "V_pt5.mat"];
colormap parula
for i = 1:4
    load(fnames(i), 'nu1', 'xu1', 'lambda_e1', 'V');
    figure(i); clf; hold on
    %saturate anything above 0.1
    lambda_e1(lambda_e1 > 0.1) = 0.1;
    contourf(xu1, nu1,lambda_e1', linspace(0,0.1,1000), 'linestyle', 'none') 
    box on
    xlim([V+0.01, 1]) %data to the left of V is meaningless, remove in this way 
    ylim([0, 10])
    c = colorbar;
    %c.Ticks = [0, 0.05,0.1];
      c.Limits = [0,0.1];
     c.TickLabels{end} = '> 0.1'; %change final tick value

end
