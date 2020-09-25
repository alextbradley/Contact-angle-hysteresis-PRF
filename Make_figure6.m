%Make figure 6 using the data file figure6data, which is generated using
%the code make_data_figure6

load('Data/figure6data.mat')
figure(1); clf; hold on
figure(2); clf; hold on
for i = 1:11
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
xlim([1e-4, 10])

figure(2);
box on
set(gca, 'XScale', 'log')
xlim([1e-4, 10])
ylim([0,hysparmax*1.2])%this cuts of the small error resulting from one sided approx and finite N