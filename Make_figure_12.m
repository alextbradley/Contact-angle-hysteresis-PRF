%Make the escape map with numerical data in Figure x in the manuscript.
%This code makes 2 plots, the first is for \lambda_{\max} < 0.1 (where
%theory works evrywhere) and for \lambda_{\max} < 0.3, in which we see
%predictions break down when corresponding equilibria are unstab;e.

%% Preliminaries
addpath('Data')
V = 0.2;
load("V_pt2.mat", 'nu2', 'xu2', 'lambda_e2', 'V', 'hysparmax_marginal','nu_marginal');
sz = 50; %size of scatter points
%% Plot 1:
figure(1); clf; hold on

%load data 
load("numerics_comparison_data_smaller_hysparmax.mat", 'hysparmax_out', 'nu_out', 'xu0_escape');

%reshape data
hysparmax_out_long = repmat(hysparmax_out, [length(nu_out),1]);
hysparmax_out_long  = hysparmax_out_long(:);
nu_out_long = repmat(nu_out,[1,length(hysparmax_out)]);
xu0_escape_long = xu0_escape';
xu0_escape_long = xu0_escape_long(:); 

%indices of always trapped or not
idx_not_always_trapped = ~isinf(xu0_escape_long);
idx_always_trapped = isinf(xu0_escape_long);
idx_always_escape = xu0_escape_long < 0; %-1 values correspond to always stuck

%Make the escape map for V = 0.3
contourf(lambda_e2, nu2,xu2', linspace(V,1,1000), 'linestyle', 'none')

%find error in predictions
prediction_error = nan(1,length(xu0_escape_long));
for i = 1:length(xu0_escape_long)
    %find nearest point in numerics
    [~,idx_nu] = min(abs(nu2 - nu_out_long(i))); 
    [~,idx_hyspar] = min(abs(lambda_e2 - hysparmax_out_long(i)));
    
    %compute error
    if ~isnan(xu2(idx_hyspar, idx_nu)) %if there is an equilibrium close
        prediction_error(i) = abs(xu0_escape_long(i) - xu2(idx_hyspar, idx_nu))/abs(xu2(idx_hyspar, idx_nu));
    end
end

%add the prediction for x_- = 0
lambda_e = linspace(0,1);
nu = 8/V^4 *(lambda_e)./(lambda_e+1).^2 .*((3*lambda_e+5)./(5*lambda_e + 5)).^4;
plot(lambda_e, nu, 'color','m', 'linewidth', 4);

%Add the numerics data
scatter(hysparmax_out_long(idx_not_always_trapped), nu_out_long(idx_not_always_trapped),...
    sz,xu0_escape_long(idx_not_always_trapped),...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2)

scatter(hysparmax_out_long(idx_always_trapped), nu_out_long(idx_always_trapped),...
    sz,'r','filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2)

scatter(hysparmax_out_long(idx_always_escape), nu_out_long(idx_always_escape),...
    sz,'g','filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2)

% %for those where error is larger than 10% and no in always trapped region (where equilibria don't exist!), use magenta outside 
% idx_big_error = (prediction_error' > 0.1) & (idx_not_always_trapped);
% scatter(hysparmax_out_long(idx_big_error), nu_out_long(idx_big_error),...
%     sz,xu0_escape_long(idx_big_error),'filled', 'MarkerEdgeColor', 'm', 'LineWidth', 2)
% 

%Plot tidy
box on
xlim([0, 0.055]) 
ylim([0, 7])

%% Plot 2:
%basically identical, with the exception of adding stability contour and
%changing axis limits
figure(2); clf; hold on

load("numerics_comparison_data_larger_hysparmax.mat", 'hysparmax_out', 'nu_out', 'xu0_escape');

%reshape data
hysparmax_out_long = repmat(hysparmax_out, [length(nu_out),1]);
hysparmax_out_long  = hysparmax_out_long(:);
nu_out_long = repmat(nu_out,[1,length(hysparmax_out)]);
xu0_escape_long = xu0_escape';
xu0_escape_long = xu0_escape_long(:); 

%indices of always trapped or not
idx_not_always_trapped = ~isinf(xu0_escape_long);
idx_always_trapped = isinf(xu0_escape_long);
idx_always_escape = xu0_escape_long < 0; %-1 values correspond to always stuck

%find errors in predictions
prediction_error = nan(1,length(xu0_escape_long));
for i = 1:length(xu0_escape_long)
    %find nearest point in numerics
    [~,idx_nu] = min(abs(nu2 - nu_out_long(i))); 
    [~,idx_hyspar] = min(abs(lambda_e2 - hysparmax_out_long(i)));
    
    %compute error
    if ~isnan(xu2(idx_hyspar, idx_nu)) %if there is an equilibrium close
        prediction_error(i) = abs(xu0_escape_long(i) - xu2(idx_hyspar, idx_nu))/abs(xu2(idx_hyspar, idx_nu));
    end
end
%Make the escape map for V = 0.3
contourf(lambda_e2, nu2,xu2', linspace(V,1,1000), 'linestyle', 'none')

%add the prediction for x_- = 0
lambda_e = linspace(0,1);
nu = 8/V^4 *(lambda_e)./(lambda_e+1).^2 .*((3*lambda_e+5)./(5*lambda_e + 5)).^4;
plot(lambda_e, nu, 'color','m', 'linewidth', 4);

%Add the numerics data
scatter(hysparmax_out_long(idx_not_always_trapped), nu_out_long(idx_not_always_trapped),...
    sz,xu0_escape_long(idx_not_always_trapped),...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2)

scatter(hysparmax_out_long(idx_always_trapped), nu_out_long(idx_always_trapped),...
    sz,'r','filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2)

%for those where error is larger than 10% and no in always trapped region (where equilibria don't exist!), use magenta outside 
idx_big_error = (prediction_error' > 0.05) & (idx_not_always_trapped);
scatter(hysparmax_out_long(idx_big_error), nu_out_long(idx_big_error),...
    sz,xu0_escape_long(idx_big_error),'filled', 'MarkerEdgeColor', 'm', 'LineWidth', 2)

%add the marginal stability contour
scatter(hysparmax_marginal, nu_marginal, 10, 'k', 'filled')

%Plot tidy
box on
xlim([0, 0.3]) 
ylim([0, 7])

