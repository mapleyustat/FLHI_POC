clear;
clc;

p0.xi = [0, 0];
p0.xo = 0;%[10];

p1.xi = [0, 1];
p1.xo = 0;%[-1];

p2.xi = [1, 0];
p2.xo = 0;%[-5];

p3.xi = [1, 1];
p3.xo = 0;%[7];

pointVector = [p0,p1,p2,p3];
kernelFunctions = [2];

[interpolant] = flhiInterpolant(pointVector, kernelFunctions);

xiVector = [0:0.01:1];

xo = zeros(length(xiVector),length(xiVector));
xi1 = zeros(length(xiVector),length(xiVector));
xi2 = zeros(length(xiVector),length(xiVector));

for xi1Index = 1:length(xiVector)
    for xi2Index = 1:length(xiVector)
        x1 = (xi1Index-1)/100;
        x2 = (xi2Index-1)/100;

        xi1(xi1Index,xi2Index) = x1;
        xi2(xi1Index,xi2Index) = x2;
        %[tmp,~] = flhiInterpolateRegion(pointVector, kernelFunctions, [x1,x2]);
        [~,tmp] = flhiInterpolate(interpolant, [x1,x2]);
        xo(xi1Index,xi2Index) = tmp(1);
    end
end

figure;
surf(xi1,xi2,xo)