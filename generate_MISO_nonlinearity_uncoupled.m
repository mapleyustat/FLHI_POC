clearvars -except FIGURE_NAME

uss1 = [-2:0.1:2];
uss2 = [-2:0.1:2];
vss1 = uss1 + 4*uss1.^2 + 1.5*uss1.^3;
vss2 = uss2 + 3*uss2.^2 + 2*uss2.^3;

fontsize = 12;

%% fig 1
h = figure('Position', [200 200 400 300]);
plot(uss1,vss1, 'k-', 'linewidth', 1.5);
ylabel('w_1', 'fontsize', fontsize);
xlabel('u_1', 'fontsize', fontsize);

%axis(gca, [min(uss) max(uss) min(vss) max(vss)]);
axis auto
set(gca, 'fontsize', fontsize);
set(gcf, 'Color', 'w');

export_fig(h, [FIGURE_NAME,'_1','.pdf'], '-q101', '-p0.01');



%% fig 2
h = figure('Position', [200 200 400 300]);
plot(uss2,vss2, 'k-', 'linewidth', 1.5);
ylabel('w_2', 'fontsize', fontsize);
xlabel('u_2', 'fontsize', fontsize);

%axis(gca, [min(uss) max(uss) min(vss) max(vss)]);
axis auto
set(gca, 'fontsize', fontsize);
set(gcf, 'Color', 'w');

export_fig(h, [FIGURE_NAME,'_2','.pdf'], '-q101', '-p0.01');