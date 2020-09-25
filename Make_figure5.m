%Make plots in figure 3 in manuscript. Data contained in
%Example_solutions_different_lambda.mat

%Note: this code does not include theta_plus as theta_plus = theta_a
%throughout

addpath('Data')

colmap = [1,0,0;
          0,1,0;
          0,0,1];
colmap = parula(4);

%Parameters
xu0 = 0.65;
xl0 = 0.45;
nu  = 4;
theta_r = 0;
tend = 10;
lambdas = [0, 0.02, 0.04];

load('Example_solutions_different_lambda', 'sols')

%initilize plots
for i = 1:5
    figure(i); clf; hold on
end

for i = 1:3
    
    %extract data from cell structure
    tout = cell2mat(sols(i,1));
    xu = cell2mat(sols(i,2));
    xl = cell2mat(sols(i,3));
    h = cell2mat(sols(i,4));
    theta_l = cell2mat(sols(i,6));
    theta_a = cell2mat(sols(i,7));
    
    if i ==2
        tmax = tout(end);
    end
    
    %plot 1: positive displacement
    figure(1);
    plot(tout, xu - xu(1), 'k', 'color', colmap(i,:), 'linewidth',3);
    plot(tout, xl - xl(1), 'k--','color', colmap(i,:), 'linewidth', 3);
    
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    xlim([1e-4, 5])
    ylim([1e-4, 1e-1])
    fig = gcf; fig.Position(3:4) = [985 211];
    box on
    xticks([1e-4, 1e-2, 1])
    
    %Plot 2 negative meniscus displacement
    figure(2);
    idx1 = ((xu - xu(1))<0); %where is xu < xu0
    idx2 = ((xl - xl(1))<0); %where is xl < xl0?
    plot(tout(idx1), xu(idx1) - xu0, 'color', colmap(i,:), 'linewidth', 3)
    plot(tout(idx2), xl(idx2) - xl0,'--','color', colmap(i,:), 'linewidth', 3)
    
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    xlim([1e-4, 5])
    ylim([-1e-2, -1e-4])
    fig = gcf; fig.Position(3:4) = [985 211];
    box on
    xticks([1e-4, 1e-2, 1])
    
    %Plot 3: pressure 
    ht = (-5*h(:, end-3) + 21*h(:,end-2) - 35*h(:,end-1) + 35*h(:, end))/16; %one sided approx to determine h at x = x_{\pm}
    hb = (-5*h(:,4) + 21*h(:,3) - 35*h(:,2) + 35*h(:,1))/16;
    p_plus = -nu*cosd(theta_a)./ht;
    p_minus = -nu*cosd(theta_l)./hb;
    
    figure(3);
    plot(tout, p_plus/nu /cosd(theta_a), 'color', colmap(i, :), 'linewidth', 3);
    plot(tout, p_minus/nu / cosd(theta_a),'--', 'color', colmap(i, :), 'linewidth', 3);
    set(gca, 'XScale', 'log')
    xlim([1e-4, 5])
    ylim([-1.1, -1])
    fig = gcf; fig.Position(3:4) = [985 211];
    box on
    xticks([1e-4, 1e-2, 1])
    
    %Plot 4: contact angle asymmetry
    figure(4)
    plot([1e-4, 5], lambdas(i)*[1,1], '--', 'color', colmap(i,:)); %dashed lines for lambda
    plot(tout, cosd(theta_l)/cosd(theta_a) -1 , 'color', colmap(i,:), 'linewidth', 3);
    set(gca, 'XScale', 'log')
    xlim([1e-4, 5])
    ylim([0, 0.05])
    fig = gcf; fig.Position(3:4) = [985 211];
    box on
    xticks([1e-4, 1e-2, 1])
    
    %Plot 5: ratio of meniscus widths
    figure(5); 
    plot(tout, hb./ht, 'color', colmap(i,:), 'linewidth', 3);
    ylim([1,1.1])
    set(gca, 'XScale', 'log')
    xlim([1e-4, 5])
    fig = gcf; fig.Position(3:4) = [985 211];
    box on
    xticks([1e-4, 1e-2, 1])

end

