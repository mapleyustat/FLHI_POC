% Organize the points in a point vector
pointVector = [];
for i = 1:length(Xo)
    for j = 1:length(Xo)
        point.xi = [i, j];
        point.xo = Xo(i,j);
        
        pointVector = [pointVector, point];
    end
end

% Create the interpolant
[interpolant] = flhiInterpolant(pointVector, INTERPOLANT_TYPE);

% Get the membership functions
MF1 = zeros(length(Xv), length(Xv));
MF2 = zeros(length(Xv), length(Xv));

for i=1:length(Xv)
    for j=1:length(Xv)
        x1 = Xv(i);
        x2 = Xv(j);
        [~, tmp] = flhiInterpolate(interpolant, [x1 x2]);
        MF1(i,j) = tmp(1);
        MF2(i,j) = tmp(3);
    end
end

% Display
[x2, x1] = meshgrid(Xv);

h = figure('Position', [200 200 400 300]);
surf(x1,x2,MF1);

%az = 120; el = 45; view(az,el);
az = 45; el = 30; view(az,el);

fontsize = 12;
xlabel('x', 'fontsize', fontsize);
ylabel('y', 'fontsize', fontsize);
zlabel('z', 'fontsize', fontsize);
set(gca, 'fontsize', fontsize);
set(gcf, 'Color', 'w');
export_fig(h, ['FLHI_MF_SURF_',FIGURE_NAME,'.pdf'], '-q101', '-p0.01');