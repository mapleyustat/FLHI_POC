% Create training data
Yt = 1./(1 + exp(-(Xt-5)));

% Organize data in a point vector
pointVector = [];
for i=1:length(Xt)
    point.xi = Xt(i);
    point.xo = Yt(i);
    
    pointVector = [pointVector, point];
end

% Create interpolant
[interpolant] = flhiInterpolant(pointVector, INTERPOLANT_TYPE);

% Interpolate validation data
Yvi = zeros(1, length(Xv));
for i=1:length(Xv)
    Yvi(i) = flhiInterpolate(interpolant, Xv(i));
end

% Interpolate inverse of training data
Xvinv = zeros(1, length(Xv));
for i=1:length(Xv)
    tmp = flhiInterpolateInverse(interpolant, Yvi(i)); % might generate multiple solutions
    Xvinv(i) = tmp(1);
end

% Display
h = figure('Position', [200 200 400 300]);
plot(Xt,Yt, 'k-', 'LineWidth', 2);
hold all;

fontsize = 12;
xlabel('x', 'fontsize', fontsize);
ylabel('y', 'fontsize', fontsize);
set(gca, 'fontsize', fontsize);
set(gcf, 'Color', 'w');
export_fig(h, ['FLHI_SISO_INV.pdf'], '-q101', '-p0.01');

h = figure('Position', [200 200 400 300]);
plot(Xv,Yvi, 'k-', 'LineWidth', 4);
hold all;
plot(Xvinv,Yvi, 'r-', 'LineWidth', 2);

fontsize = 12;
xlabel('x', 'fontsize', fontsize);
ylabel('y', 'fontsize', fontsize);
set(gca, 'fontsize', fontsize);
set(gcf, 'Color', 'w');
export_fig(h, ['FLHI_SISO_INV_INTERPOLATED_',FIGURE_NAME,'.pdf'], '-q101', '-p0.01');

% Error analysis
residual = Xv - Xvinv;
disp(['RESIDUAL FLHI SISO INVERSE - ', FIGURE_NAME]);
disp('mean - std - min|| - max||');
[mean(residual) std(residual) min(abs(residual)) max(abs(residual))]