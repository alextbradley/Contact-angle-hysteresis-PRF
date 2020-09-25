%Produce equilibria maps in (\lambda_e, \nu) space with colour according to
%X_+for V = 0.1,0.2, 0.3, 0.4, 0.5 
figpref(4);
addpath('./Data')
fnames = ["V_pt1.mat","V_pt2.mat","V_pt3.mat", "V_pt4.mat", "V_pt5.mat"];
colormap parula
for i = 1:5
    load(fnames(i), 'nu2', 'xu2', 'lambda_e2', 'V');
    figure(i); clf; hold on
    contourf(lambda_e2, nu2,xu2', linspace(V,1,1000), 'linestyle', 'none') 
    box on
    xlim([0, 0.1]) 
    ylim([0, 10])
    
    %colorbar
%     c = colorbar; %tick labels might appear different from paper
%     c.Limits = [V,1];

    %add the prediction for x_- = 0
    lambda_e = linspace(0,1);
    nu = 8/V^4 *(lambda_e)./(lambda_e+1).^2 .*((3*lambda_e+5)./(5*lambda_e + 5)).^4;
    hold on
    plot(lambda_e, nu, 'color','m', 'linewidth', 4);
    
    
end