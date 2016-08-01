% Create interpolation training data
x1 = zeros(length(Xt),length(Yt));
x2 = zeros(length(Xt),length(Yt));
x3 = zeros(length(Xt),length(Yt));
x4 = zeros(length(Xt),length(Yt));
pointVector = [];
for i=1:length(Xt)
    for j=1:length(Yt)
        x1(i,j) = Xt(i);
        x2(i,j) = Yt(j);
        
        x = x1(i,j);
        y = x2(i,j);
        
        % Binh and Korn function
        x3(i,j) = 4*x^2 + 4*y^2;
        x4(i,j) = (x-5)^2 + (y - 5)^2;
        
        % Organize data in a point vetor
        point.xi = [x1(i,j), x2(i,j)];
        point.xo = [x3(i,j), x4(i,j)];
        pointVector = [pointVector, point];
    end
end

% Create the interpolant
[interpolant] = flhiInterpolant(pointVector, INTERPOLANT_TYPE);

% Interpolate validation data
xi1 = zeros(length(Xv),length(Yv));
xi2 = zeros(length(Xv),length(Yv));
xi3 = zeros(length(Xv),length(Yv));
xi4 = zeros(length(Xv),length(Yv));
xr3 = zeros(length(Xv),length(Yv));
xr4 = zeros(length(Xv),length(Yv));

for i = 1:length(Xv)
    for j = 1:length(Yv)
        xi1(i,j) = Xv(i);
        xi2(i,j) = Yv(j);
        [tmp,~] = flhiInterpolate(interpolant, [xi1(i,j),xi2(i,j)]);
        xi3(i,j) = tmp(1);
        xi4(i,j) = tmp(2);
        
        x = xi1(i,j);
        y = xi2(i,j);
        
        % Binh and Korn function
        xr3(i,j) = 4*x^2 + 4*y^2;
        xr4(i,j) = (x-5)^2 + (y - 5)^2;
    end
end

% Display
h = figure('Position', [200 200 450 350]);
surf(x1,x2,x3);
hold all;
surf(x1,x2,x4);

fontsize = 12;
xlabel('x', 'fontsize', fontsize);
ylabel('y', 'fontsize', fontsize);
zlabel('z', 'fontsize', fontsize);
set(gca, 'fontsize', fontsize);
set(gcf, 'Color', 'w');
export_fig(h, ['FLHI_MIMO_DATA.pdf'], '-q101', '-p0.01');

%
h = figure('Position', [200 200 450 350]);
surf(xi1,xi2,xi3);
hold all;
surf(xi1,xi2,xi4);

fontsize = 12;
xlabel('x', 'fontsize', fontsize);
ylabel('y', 'fontsize', fontsize);
zlabel('z', 'fontsize', fontsize);
set(gca, 'fontsize', fontsize);
set(gcf, 'Color', 'w');
export_fig(h, ['FLHI_MIMO_INTERPOLATED_',FIGURE_NAME,'.pdf'], '-q101', '-p0.01');

% Error analysis
residual1 = xr3 - xi3; residual1 = residual1(:);
residual2 = xr4 - xi4; residual2 = residual2(:);
disp(['RESIDUAL FLHI MIMO - ', FIGURE_NAME]);
disp('mean - std - min - max');
[mean(residual1) std(residual1) min(abs(residual1)) max(abs(residual1))]
[mean(residual2) std(residual2) min(abs(residual2)) max(abs(residual2))]
