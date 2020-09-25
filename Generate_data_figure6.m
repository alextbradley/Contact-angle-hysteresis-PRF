%Generate data for figure 6 in manuscript. Note that this code also
%produces the figure. Warning! This data takes O(1 hour) to generate!

%intialize figures
figure(1);clf; hold on
figure(2); clf; hold on

%Parameters
N = 50; %large numbers of grid points not recommended
hysparmax = 0.02;
nu = 2;
V = 0.2;
tend = 10;
logout = 1;
theta_r = 0;
theta_a = acosd(cosd(theta_r)/(1 + hysparmax));
xu0s = [0.5, 0.55, 0.6, 0.65, 0.7, 0.725, 0.75,0.77, 0.8, 0.85, 0.9];
data = cell(length(xu0s),5);
for k = 1:length(xu0s)
    xu0 = xu0s(k);

    [v, t, t_events_functions, ispinned] = run_droplet(nu,hysparmax, theta_r, V, xu0,N, tend, logout);
   
    %Process the solution
    sz = size(v);
    Nruns = size(ispinned);
    Nruns = Nruns(1); %number of runs in the integration
    Nstepsout = round(sz(1)/Nruns); %number of timesteps per integration
    
    xl = v(:,end-1);
    xu = v(:,end);
    w = xu - xl;
    dz = 1/N;
    pxl = (7*v(:,1) - 33*v(:,2) + 62*v(:,3) - 58*v(:,4) + 27*v(:,5) - 5*v(:,6))/2 /dz^4 ./w.^5;%pressure using one sided derivative
    uxl = (35*v(:,1) - 35*v(:,2) + 21*v(:,3) - 5*v(:,4))/16 ; %u at Z = 0
    pxu = (7*v(:,N) - 33*v(:,N-1) + 62*v(:,N-2) - 58*v(:,N-3) + 27*v(:,N-4) - 5*v(:,N-5))/2 /dz^4 ./w.^5;%pressure using one sided derivative
    uxu = (35*v(:,N) - 35*v(:,N-1) + 21*v(:,N-2) - 5*v(:,N-3))/16 ; %u at Z = 1
    
    theta_l = acosd(-pxl.*uxl.*cosd(theta_a)/nu ./w); %this should work at all times, but loses accuracy because of low N and use of central FD in numerics. So only use this formulation when pinned and we don't know explicit value of contact angles
    theta_u = acosd(-pxu.*uxu.*cosd(theta_a)/nu ./w);
    for i = 1:Nruns
        switch ispinned(i,1) %cases for pinning conditions at xl
            case 1 %receding angle
                theta_l((i-1)*Nstepsout + 1: i*Nstepsout) =acosd(cosd(theta_r) + 1e-3);
            case 2
                theta_l((i-1)*Nstepsout + 1: i*Nstepsout) = theta_a;
        end
        
        switch ispinned(i,2) %cases for pinning conditions at xu
            
            case 1 %receding angle
                theta_u((i-1)*Nstepsout + 1: i*Nstepsout) = theta_r;
            case 2
                theta_u((i-1)*Nstepsout + 1: i*Nstepsout) = theta_a;
        end
    end
    
    %add to data file
    data{k,1} = xu0;
    data{k,2} = t;
    data{k,3} = xu;
    data{k,4} = cosd(theta_l)./cosd(theta_u) - 1;
    data{k,5} = ispinned;


end

%save data 
save('Data/figure6data.mat', 'data', 'hysparmax', 'nu', 'theta_r', 'theta_a')
%% plot

for i = 1:length(xu0s)
    t = cell2mat(data(i,2));
    xu = cell2mat(data(i,3));
    lambda = cell2mat(data(i,4));
    ispinned =cell2mat(data(i,5));
    figure(1);
    if ispinned(end,1) == 1
        plot(t, xu, 'color', [209, 8,8]/255, 'linewidth', 3)
    else
        plot(t, xu, 'color', [253, 191,121]/255, 'linewidth', 3)
    end
    
    
    figure(2);
    if ispinned(end,1) == 1
        plot(t, lambda,'color', [209, 8,8]/255, 'linewidth', 3)
    else
        plot(t, lambda,'color', [253, 191,121]/255, 'linewidth', 3)
    end
end

figure(1);
box on
set(gca, 'XScale', 'log')
xlim([1e-4, tend])

figure(2);
box on
set(gca, 'XScale', 'log')
xlim([1e-4, tend])
ylim([0,hysparmax*1.2])%this cuts of the small error resulting from one sided approx and finite N