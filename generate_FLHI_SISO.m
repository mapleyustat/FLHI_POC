% Create training data
%Xt = -3:0.1:3;
Zt = peaks(Xt,0);

% Create validation data
%Xv = -3:0.01:3;
%Xv = 2.95:0.05:3;
Zv = peaks(Xv,0);

% Organize training data in a point vector
pointVector = [];
for i=1:length(Xt)
    point.xi = Xt(i);
    point.xo = Zt(i);
    
    pointVector = [pointVector, point];
end

% Create interpolant
[interpolant] = flhiInterpolant(pointVector, INTERPOLANT_TYPE);

% Interpolate validation data
Zvi = zeros(1, length(Zv));
for i=1:length(Xv)
    Zvi(i) = flhiInterpolate(interpolant, Xv(i));
end

% Display
h = figure('Position', [200 200 400 300]);
plot(Xt,Zt, 'k--', 'LineWidth', 1);
hold all;
plot(Xt,Zt, 'ko', 'LineWidth', 2);

fontsize = 12;
xlabel('x', 'fontsize', fontsize);
ylabel('z', 'fontsize', fontsize);
set(gca, 'fontsize', fontsize);
set(gcf, 'Color', 'w');
export_fig(h, ['FLHI_SISO_DATA.pdf'], '-q101', '-p0.01');

h = figure('Position', [200 200 400 300]);
plot(Xv,Zv, 'k', 'LineWidth', 2);
hold all;
plot(Xv,Zvi, 'b', 'LineWidth', 2);

legend('Real', 'FLHI', 'Location', 'southeast');
fontsize = 12;
xlabel('x', 'fontsize', fontsize);
ylabel('z', 'fontsize', fontsize);
set(gca, 'fontsize', fontsize);
set(gcf, 'Color', 'w');
export_fig(h, ['FLHI_SISO_INTERPOLATED_',FIGURE_NAME,'.pdf'], '-q101', '-p0.01');

% Error analysis
residual = Zv - Zvi;
disp(['RESIDUAL FLHI SISO - ', FIGURE_NAME]);
disp('mean - std - min|| - max||');
[mean(residual) std(residual) min(abs(residual)) max(abs(residual))]