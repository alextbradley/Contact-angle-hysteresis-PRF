%Make figure 3 in manuscript.
clear
colors = [80, 62, 182; 
        65, 145, 255;
        0, 200, 197;
        255, 174, 103]/255;
figpref(4);
%% (a) Petrov and Petrov 1991 data on top of our model
figure(1); clf;

%load everything from petrov petrov data file
load('data/PetrovPetrov1991Fig6Data.mat');

%add our function as base
subplot(1,2,1); hold on
xdot = linspace(-30,30); %units of microns/sec
theta = theta_a*(xdot > 0) + theta_r*(xdot < 0);
plot(xdot, theta_a*ones(size(xdot)), 'k--')
plot(xdot, theta_r*ones(size(xdot)), 'k--')
plot(xdot, theta, 'linewidth', 3, 'color', colors(1,:))


%add PetrovPetrov data
%(scale of 10 because Petrov data recording in 10^-3 cm/s)
plot(v_advancing*10, theta_advancing_dynamic,'o', 'color', colors(4,:))
plot(-v_receding*10, theta_receding_dynamic,'o', 'color',  colors(4,:))

%tidy up
box on
xlim([-30, 30]);
xlabel('$\dot{x}_+$~ ($\mu \mathrm{m~s}^{-1}$)', 'interpreter', 'latex') 
ylabel('$\theta_{a,r}$ (degrees)', 'interpreter', 'latex')
txta = text(-40, 80, '(a)', 'interpreter', 'latex', 'FontSize', 14);
txtadv = text(-20, 68.5, '$\theta_a = 70.5{}^\circ$', 'FontSize', 12, 'interpreter', 'latex');
txtrec = text(10, 40.5, '$\theta_r = 38.5{}^\circ$', 'FontSize', 12, 'interpreter', 'latex');

%% (b) lambda_max as a function of theta_a - theta_r 
%for different theta_r
theta_rs = [45, 30, 15, 0];
subplot(1,2,2); hold on
for i = 1:4
    theta_r = theta_rs(i);
    theta_a = linspace(theta_r, 70);
    lambda_max = cosd(theta_r) ./ cosd(theta_a) -1;
    plot(theta_a - theta_r, lambda_max, 'linewidth', 3, 'color', colors(i,:));
end

box on
xlim([0, 50]);
ylim([0, 0.5]);
xlabel('$\theta_a - \theta_r$ (degrees)', 'interpreter', 'latex')
ylabel('$\lambda_{\max}$', 'interpreter', 'latex')
txtb = text(-9, 0.5, '(b)', 'interpreter', 'latex', 'FontSize', 14);
txt = text(8, 0.47, '$\theta_r =$', 'interpreter', 'latex', 'FontSize', 12);
txt45 = text(13, 0.47, '$45^\circ$', 'interpreter', 'latex', 'FontSize', 12, 'color', colors(1,:));
txt30 = text(21, 0.47, '$30^\circ$', 'interpreter', 'latex', 'FontSize', 12, 'color', colors(2,:));
txt15 = text(31, 0.47, '$15^\circ$', 'interpreter', 'latex', 'FontSize', 12, 'color', colors(3,:));
txt0 = text(45, 0.47, '$0^\circ$', 'interpreter', 'latex', 'FontSize', 12, 'color', colors(4,:));

fig = gcf; 
fig.Position(3:4) = [1250, 340];
shg