%Make figure 3 in manuscript.
clear
colors = [80, 62, 182; 
        65, 145, 255;
        0, 200, 197;
        255, 174, 103]/255;
%% (a) Petrov and Petrov 1991 data on top of our model
figure(1); clf;

%load everything from petrov petrov data file
load('data/PetrovPetrov1991Fig6Data.mat');

%add our function as base
subplot(1,2,1); hold on
xdot = linspace(-30,30); %units of microns/sec
theta = theta_a*(xdot > 0) + theta_r*(xdot < 0);
plot(xdot, theta, 'linewidth', 3, 'color', colors(1,:))

%add PetrovPetrov data
%(scale of 10 because Petrov data recording in 10^-3 cm/s)
plot(xdot, theta_a*ones(size(xdot)), 'k--')
plot(xdot, theta_r*ones(size(xdot)), 'k--')
plot(v_advancing*10, theta_advancing_dynamic,'o', 'color', colors(4,:))
plot(-v_receding*10, theta_receding_dynamic,'o', 'color',  colors(4,:))

%tidy up
box on
xlim([-30, 30]);
xlabel('$\dot{x}_2$~ ($\mu \mathrm{m~s}^{-1}$)', 'interpreter', 'latex') 
ylabel('$\theta_{a,r}$ (${}^\circ$)', 'interpreter', 'latex')
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
xlabel('$\theta_a - \theta_r$ (${}^\circ$)', 'interpreter', 'latex')
ylabel('$\lambda_{\max}$', 'interpreter', 'latex')