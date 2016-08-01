function [error] = fobjRegionInverseSearch(xi, region, kernelFunctions, xoExpected)
%Objective function to search which input coordinates lead to the expected
%output coordinates in a region
%   Something, something... dark side.

xo = flhiInterpolateRegion(region.pointVector, kernelFunctions, xi);

error = sum((xo - xoExpected).^2);

end

