%Produce regime diagram for hyspar_e = 0.3 and 0.05
%% Preliminaries
clear; close all
addpath('./Data')
addpath('Auxillary_functions')
figpref(4);
%% Load data and plot
count = 1;
dTs = [0.05, 0.3];
for fname = ["lambda_e_pt05.mat", "lambda_e_pt3.mat"]
    figure(count); clf; hold on
    load(fname, 'X', 'Y')

    pl = fill(X,Y, 'r', 'edgecolor', 'r');
    
    % Add x_- = 0 theory
    dT = dTs(count);
    V_int = 4*dT.*(3*dT + 5)./5./(dT+1) ./(4*dT+3); %Volume where cross over to regime II
    vv = linspace(V_int,1);
    nu = 8./(dT+1) ./vv.^4 .*(1-1./(dT+1)) .*(3/5 + 2/5./(dT+1)).^4;
    plot(vv, nu,'k--', 'linewidth', 3);
    ylim([0, 100])
    box on
    xlabel('$V$', 'interpreter', 'latex')
    ylabel('$\nu$', 'interpreter', 'latex')
    
    %make insets
    ax = axes('Position',[.6 .6 .28 .28]);
    ax.FontSize = 16;
    pl = fill(X,Y, 'r', 'edgecolor', 'r');
    ylim([0, 10])
    box on
    count = count + 1;
    
    
end

%% add the insets
figure(1);
pars = [40, 0.01; %medium nu, close to base
        80, 0.01; %large nu, close to base
        2, 0.1; % large volume, small nu
        20,0.8]; %low volume, close to end

axloc = [.3, .4;
         .3, .7;
         .4, .15;
         .14, .22];
     
axsz = [.15, .07]; %size of the axes

for i = 1:4
    ax = axes('Position',[axloc(i,:),axsz]);
    pleq_schematic(1.05, pars(i,1), pars(i,2)) %first argument is lambda + 1
end
    
%figure 2:
figure(2);
dT = 0.3;
pars = [6, 0.7; %low volume, near tips
        50, 0.25; %close to contact midway up
        80, 0.2; %close to contact halfway up
        5, 0.1; %high volume, close to base
        40, 0.1; %midway up close to base
        80, 0.1]; %midway up, close to base
    
 axloc = [.15,.15;
          .15, .5;
          .15, .8;
          .52, .15;
          .35, .5;
          .35, .8];
      
for i = 1:6
    ax = axes('Position',[axloc(i,:),axsz]);
    pleq_schematic(1.3, pars(i,1), pars(i,2))  %first argument is lambda + 1
end
          
          
        
        
        

