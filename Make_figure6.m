%Make figure 6 using the data file figure6data, which is generated using
%the code make_data_figure6

load('Data/figure6data.mat')
figure(1); clf; hold on

%make colormaps
colmap_trapped = zeros(7,3);
colmap_trapped(:,1) = linspace(0.5,1,7);

colmap_escape = zeros(4,3);
colmap_escape(:,end) = flip(linspace(0.3,1,4));
colmap = [colmap_trapped; colmap_escape];

colmap_full = parula(18);
colmap = [colmap_full(1:7,:); colmap_full(end-3:end, :)];


for i = 1:11
    t = cell2mat(data(i,2));
    xu = cell2mat(data(i,3));
    lambda = cell2mat(data(i,4));
    ispinned =cell2mat(data(i,5));
    subplot(1,2,1); hold on
    if ispinned(end,1) == 1
        plot(t, xu, 'color', colmap(i,:), 'linewidth', 3)
    else
        plot(t, xu, 'color', colmap(i,:), 'linewidth', 3)
    end
    
    
    subplot(1,2,2);  hold on
    if ispinned(end,1) == 1
        plot(t, lambda,'color', colmap(i,:), 'linewidth', 3)
    else
        plot(t, lambda,'color', colmap(i,:), 'linewidth', 3)
    end
end

subplot(1,2,1);
box on
set(gca, 'XScale', 'log')
xlim([1e-4, 10])
xlabel('$t$', 'interpreter', 'latex');
ylabel('$x_2(t)$', 'interpreter', 'latex');

subplot(1,2,2);
box on
set(gca, 'XScale', 'log')
xlim([1e-5, 10])
ylim([0,hysparmax*1.2])%this cuts of the small error resulting from one sided approx and finite N
xlabel('$t$', 'interpreter', 'latex');
ylabel('$\lambda(t)$', 'interpreter', 'latex');

fig = gcf;
fig.Position(3:4) = [1200, 420];