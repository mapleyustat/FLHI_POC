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

xiVector = [minx:maxx];

pointVector = [];

for i = 1:length(xiVector)
    for j = 1:length(xiVector);
        point.xi = [xiVector(i),xiVector(j)];
        point.xo = peaks(xiVector(i),xiVector(j));
        
        pointVector = [pointVector, point];
    end
end

% Interpolate
[interpolant] = flhiInterpolant(pointVector, 2);

xiVector = [minx:0.1:maxx];

xo = zeros(length(xiVector),length(xiVector));
xi1 = zeros(length(xiVector),length(xiVector));
xi2 = zeros(length(xiVector),length(xiVector));

for i = 1:length(xiVector)
    for j = 1:length(xiVector)
        xi1(i,j) = xiVector(i);
        xi2(i,j) = xiVector(j);
        [tmp,~] = flhiInterpolate(interpolant, [xi1(i,j),xi2(i,j)]);
        %[~,tmp] = flhiInterpolate(interpolant, [xi1(i,j),xi2(i,j)]);
        xo(i,j) = sum(tmp);
    end
end

figure;
surf(xi1,xi2,xo);

% figure;
% I = (xi2 >= -3 & xi2 < -2) & (xi1 >= 2 & xi1 < 3);
% xii1 = reshape(xi1(find(I)), max(sum(I)), max(sum(I')));
% xii2 = reshape(xi2(find(I)), max(sum(I)), max(sum(I')));
% xio = reshape(xo(find(I)), max(sum(I)), max(sum(I')));
% surf(xii1, xii2, xio);

xi = flhiInterpolateInverse(interpolant, 5.859)

