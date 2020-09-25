%Make solvability contour plots:
%Plot 1: 6 subplots mapping equilibria in (nu,V) space for nu < 10. Each
%subplot corresponds to a different value of \lambdae as follows: 0.01,
%0.05, 0.1, 0.2, 0.3, 0.4. The black line indicates the zero contour of the
%solvability condition S (defined in manuscript) and the coloured vertical
%line indicates the line along which the growth rate sigma is computed (V =
%0.3 in each case), as is plotted in plot 2.
%
%Plot 2: plot the growth rate of the perturbation as a function of \nu,
%taken along the coloured line in each of the subplots of (a). Growth rates
%are determined by solving the BVP, but for $\lambdae = 0.05, we also plot
%the growth rate obtained by fitting to the data at early times.

%NOTE 1: this code applies a very crude search mechanism to get initial
%guesses for growth rate. The values of initguess and increment associated
%with this search are obtained by trial and error(!) and it is not
%recommended that they are changed (there is no safety mechanism to prevent
%search continuing ad infinitum)

%NOTE 2: This code take O(minutes) to run

addpath('Data')
%% Plot 1
%load data
load('Small_nu_solvability_values.mat')

%   %define colormaps
xx = linspace(0,1);
cmap1 = zeros(100,3); %fade to red
cmap1(:,1) = xx.^(1/2); %change exponent to control shading.
cmap1(:,1) = 1; %red everywhere
figure(1); clf;
for i = 1:6 %don't include \lmabdae = 0.5
    subplot(2,3,i)
    S = cell2mat(dataR1(1,i));
    V_ = cell2mat(dataR1(3,i));
    nu_ = cell2mat(dataR1(2,i));
    [c,h] = contourf(V_, nu_, log(abs(S')), -5:0.01:5);
    h.LineStyle = 'none';
    hold on 
    contour(V_, nu_,abs(S)', [1e-2,1e-2],'k'); %add the zero contour (choose a level close to zero). We touch up this plot in a graphics editor to join these contours. 
    colormap(cmap1);
    %colorbar
    
    xlim([0.02,1])
    ylim([0 10]);
end
tic
%% Plot 2
%first run 0.05 case. uses numerics to guide initial guesses (unlike other
%cases below)
V = 0.3;
hyspar = 0.05;
nus = linspace(1,10,40);

%initialize storage
nus_with_equilibrium = nan(6,length(nus)); 
sigma_numerics = nan(6,length(nus));
sigma_BVP = nan(6,length(nus));

for i = 1:length(nus)
    sigma_numerics(2,i) = growthrate_numerics(V, hyspar, nus(i), 40); %final argument number of grid points in FD scheme
    initguess = sigma_numerics(2,i);
    if ~isnan(initguess)
        nus_with_equilibrium(2,i) = nus(i); %keep track of those nu values with an equilibrium associated
        sigma_BVP(2,i)      = growthrate_BVP(V, hyspar, nus(i),initguess);
    end

end

subplot(2,3,2); hold on
plot(V*ones(1, length(nus_with_equilibrium(2,:))), nus_with_equilibrium(2,:), 'k--', 'linewidth', 3)

% figure(2); clf;
% plot(nus_with_equilibrium(2,:), sigma_numerics(2,:), 'ko'); hold on
% plot(nus_with_equilibrium(2,:), sigma_BVP(2,:), 'k--');
fprintf('finished lambdae = 0.05 \n');
toc
%% Other values of lambdae
init_guesses = [-5000,nan, -150, -25, -5, 5]; %second entry corresponds to lambdae = 0.05 (already completed above)
increments = [100,nan, 10, 0.2, 0.2, 0.2];
lambdeS = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4];

for k = [1,3,4,5,6] %don't include lambda = 0.05 in your loop
    lambda = lambdeS(k);
    for i = 1:length(nus)
        initguess = init_guesses(k);
        sigma_BVP(k,i) = growthrate_BVP(V, lambda, nus(i),initguess);
        
        %this value of nu corresponds to an equilibrium if you returned not nan
        if ~isnan(sigma_BVP(k,i))
            nus_with_equilibrium(k,i) = nus(i);
        end
        
        % if you return a zero, try a different initial guess
        while abs(sigma_BVP(k,i)) < 1e-3 %hopefully we don't hit any zeros actually very close to zero!
            initguess = initguess + increments(k);
            sigma_BVP(k,i) = growthrate_BVP(V, lambda, nus(i),initguess);
        end
        
    end
    toc
    fprintf('fininshed lambdae = %.2f \n', lambda)
end


%% Plots:
cmap = parula(length(lambdeS) + 1);

figure(1);
%add dashed lines on figure 1 
for k = 1:length(lambdeS)
    subplot(2,3,k); hold on
    plot(V*ones(1, length(nus_with_equilibrium(k,:))), nus_with_equilibrium(k,:),...
        '--', 'linewidth', 3, 'color', cmap(k,:))
    xlim([0,1])
    xticks([0, 0.5, 1])
end
fig = gcf; fig.Position(3:4) = [878 545];
%%
figure(2); clf
hold on
%add grey dashed line for sigma = 0
plot([0,10], [0,0],'--', 'color', [0.5, 0.5, 0.5])
for k = 4:length(lambdeS)
    plot(nus_with_equilibrium(k,:), sigma_BVP(k,:), 'linewidth', 3, 'color', cmap(k,:))
end
ylim([-20,10])
xlim([0,10])
box on
fig = gcf; fig.Position(3:4) = [500 228];

% repeat with more outzoom
figure(3); clf
hold on
%add grey dashed line for sigma = 0

for k = 1:3
    plot(nus_with_equilibrium(k,:), sigma_BVP(k,:), 'linewidth', 3, 'color', cmap(k,:))
end
xlim([0,10])
box on
fig = gcf; fig.Position(3:4) = [500 228];
ylim([-500,0])
yticks([-500,-250,0])
xlim([0,10])
%add numerics for lambda = 0.05
plot(nus_with_equilibrium(2,:), sigma_numerics(2,:),'--', 'linewidth', 3, 'color', cmap(2,:))

%add inset with lambde = 0.01 on
axes('Position',[.23 .23 .2 .3])
box on
plot(nus_with_equilibrium(1,:), sigma_BVP(1,:), 'linewidth', 3, 'color', cmap(1,:))
xlim([0,10])
ylim([-6000, 0])
