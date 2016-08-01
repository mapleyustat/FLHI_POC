clearvars -except FIGURE_NAME

uss1 = [-2:0.1:2];
uss2 = [-2:0.1:2];
%vss1 = uss1 + 4*uss1.^2 + 1.5*uss1.^3;
%vss2 = uss2 + 3*uss2.^2 + 2*uss2.^3;

[USS1,USS2] = meshgrid(uss1,uss2);
VSS1 = USS1 + 4*USS1.^2 + 1.5*USS1.^3;
VSS2 = USS2 + 3*USS2.^2 + 2*USS2.^3;

fontsize = 12;

%% fig 1
h = figure('Position', [200 200 400 300]);
s = surf(USS1,USS2,VSS1);

set(s,'facecolor','none')
%colormap(gray)

zlabel('w_1', 'fontsize', fontsize);
ylabel('u_2', 'fontsize', fontsize);
xlabel('u_1', 'fontsize', fontsize);

axis auto
set(gca, 'fontsize', fontsize);
set(gcf, 'Color', 'w');

export_fig(h, [FIGURE_NAME,'_1','.pdf'], '-q101', '-p0.01');

%% fig 2
h = figure('Position', [200 200 400 300]);
s = surf(USS1,USS2,VSS2);

set(s,'facecolor','none')
%colormap(gray)

zlabel('w_2', 'fontsize', fontsize);
ylabel('u_2', 'fontsize', fontsize);
xlabel('u_1', 'fontsize', fontsize);

axis auto
set(gca, 'fontsize', fontsize);
set(gcf, 'Color', 'w');

export_fig(h, [FIGURE_NAME,'_2','.pdf'], '-q101', '-p0.01');