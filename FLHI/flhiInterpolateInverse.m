function [xi, fval] = flhiInterpolateInverse(interpolant, xo)

ximax = 1;
ximin = 0;

ub = ones(1, interpolant.inputDimensions)*ximax;
lb = ones(1, interpolant.inputDimensions)*ximin;
xiInitial = ones(1, interpolant.inputDimensions)*((ximax-ximin)/2);

optoptions = optimset('Display','off','Algorithm','interior-point','TolFun',1e-10);

fval = zeros(interpolant.regionCount, 1);
xi = zeros(interpolant.regionCount, interpolant.inputDimensions);
xiIndex = 1;

% TODO: make the region search more efficient - maybe kd-tree based?
for regionIndex = 1:interpolant.regionCount
    region = interpolant.regionVector(regionIndex);
    
    isInUpperBound = all(xo <= region.xoMaximum);
    isInLowerBound = all(xo >= region.xoMinimum);
    isInRegion = isInUpperBound && isInLowerBound;
    
    if isInRegion
         [xiFound,fval(xiIndex)] = fmincon(@(xi)fobjRegionInverseSearch(xi,region,interpolant.kernelFunctions,xo),xiInitial,[],[],[],[],lb,ub,[],optoptions);
         
         xi(xiIndex, :) = xiFound.*interpolant.xiStep + region.basePoint.xi;
         
         xiIndex = xiIndex + 1;
    end
end

xi(xiIndex:end, :) = [];
fval(xiIndex:end, :) = [];

end

