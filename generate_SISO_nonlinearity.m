clearvars -except FIGURE_NAME

uss = [-0.15:0.01:0.15];
vss = 1.04.*uss - 14.11.*uss.^2 - 16.72.*uss.^3 + 562.75.*uss.^4;

fontsize = 12;
h = figure('Position', [200 200 400 300]);
plot(uss,vss, 'k-', 'linewidth', 1.5);
ylabel('w', 'fontsize', fontsize);
xlabel('u', 'fontsize', fontsize);

%axis(gca, [min(uss) max(uss) min(vss) max(vss)]);
axis auto
set(gca, 'fontsize', fontsize);
set(gcf, 'Color', 'w');

export_fig(h, [FIGURE_NAME,'.pdf'], '-q101', '-p0.01');
