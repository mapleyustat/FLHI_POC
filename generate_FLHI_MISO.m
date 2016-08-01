% Create interpolation training data
[x1,x2] = meshgrid(Xt);
x3 = peaks(x1,x2);

% Organize data in a point vetor
pointVector = [];
for i = 1:length(Xt)
    for j = 1:length(Xt);
        point.xi = [Xt(i),Xt(j)];
        point.xo = peaks(Xt(i),Xt(j));
        
        pointVector = [pointVector, point];
    end
end

% Create the interpolant
[interpolant] = flhiInterpolant(pointVector, INTERPOLANT_TYPE);

% Interpolate validation data
xo = zeros(length(Xv),length(Xv));
xi1 = zeros(length(Xv),length(Xv));
xi2 = zeros(length(Xv),length(Xv));

xoo = zeros(length(Xv),length(Xv));

for i = 1:length(Xv)
    for j = 1:length(Xv)
        xi1(i,j) = Xv(i);
        xi2(i,j) = Xv(j);
        [tmp,~] = flhiInterpolate(interpolant, [xi1(i,j),xi2(i,j)]);
        %[~,tmp] = flhiInterpolate(interpolant, [xi1(i,j),xi2(i,j)]);
        xo(i,j) = tmp(1);
        xoo(i,j) = peaks(xi1(i,j),xi2(i,j));
    end
end

% Display
h = figure('Position', [200 200 450 350]);
surf(x1,x2,x3);

fontsize = 12;
xlabel('x', 'fontsize', fontsize);
ylabel('y', 'fontsize', fontsize);
zlabel('z', 'fontsize', fontsize);
set(gca, 'fontsize', fontsize);
set(gcf, 'Color', 'w');
export_fig(h, ['FLHI_MISO_DATA.pdf'], '-q101', '-p0.01');

%
h = figure('Position', [200 200 450 350]);
surf(xi1,xi2,xoo);

fontsize = 12;
xlabel('x', 'fontsize', fontsize);
ylabel('y', 'fontsize', fontsize);
zlabel('z', 'fontsize', fontsize);
set(gca, 'fontsize', fontsize);
set(gcf, 'Color', 'w');
export_fig(h, ['FLHI_MISO_REAL.pdf'], '-q101', '-p0.01');

%
h = figure('Position', [200 200 450 350]);
surf(xi1,xi2,xo);

fontsize = 12;
xlabel('x', 'fontsize', fontsize);
ylabel('y', 'fontsize', fontsize);
zlabel('z', 'fontsize', fontsize);
set(gca, 'fontsize', fontsize);
set(gcf, 'Color', 'w');
export_fig(h, ['FLHI_MISO_INTERPOLATED_',FIGURE_NAME,'.pdf'], '-q101', '-p0.01');

% Error analysis
residual = xoo - xo;
residual = residual(:);
disp(['RESIDUAL FLHI MISO - ', FIGURE_NAME]);
disp('mean - std - min - max');
[mean(residual) std(residual) min(abs(residual)) max(abs(residual))]

