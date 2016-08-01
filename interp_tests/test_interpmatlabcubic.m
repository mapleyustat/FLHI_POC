clear;
clc;

% Create interpolation data
minx = -3;
maxx = 3;
[x1,x2] = meshgrid(minx:1:maxx);
x3 = peaks(x1,x2);

figure;
surf(x1,x2,x3);
%plot3(x1,x2,x3,'or');

% Interpolate
[Xq,Yq] = meshgrid(-3:0.1:3);

Vq = interp2(x1,x2,x3,Xq,Yq,'spline');

figure;
surf(Xq,Yq,Vq);
return

% figure;
% I = (xi2 >= -3 & xi2 < -2) & (xi1 >= 2 & xi1 < 3);
% xii1 = reshape(xi1(find(I)), max(sum(I)), max(sum(I')));
% xii2 = reshape(xi2(find(I)), max(sum(I)), max(sum(I')));
% xio = reshape(xo(find(I)), max(sum(I)), max(sum(I')));
% surf(xii1, xii2, xio);

xi = flhiInterpolateInverse(interpolant, 5.859)

