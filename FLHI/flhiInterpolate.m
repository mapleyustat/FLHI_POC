function [xo, membershipVector] = flhiInterpolate(interpolant, xi)
%Interpolates a position 'xi' given a data set 'interpolant'
%   Something, something... dark side.

% TODO: Verify if inputDimensions from xi is equal to interpolant.inputDimensions
xiIndex = zeros(1, interpolant.inputDimensions);
xiHipercube = zeros(1, interpolant.inputDimensions);

% translate coordinates to logical indices for table look up
for i = 1:interpolant.inputDimensions
    % TODO: BOUNDS CHECK!
    xiIndex(i) = floor((xi(i) - interpolant.xiMinimum(i))/interpolant.xiStep(i));
    %xiIndex(i) = floor((xi(i) - interpolant.xiMinimum(i))/interpolant.xiStep(i) + 1000*eps); % TODO: more reliable against numerical errors?
end

regionIndex = 1; % (start at 1, because Matlab...)
len = 1;
for i = 1:interpolant.inputDimensions
    % if a point is located at the upper edge, make sure an existing region
    % will be mapped
    if xi(i) == interpolant.xiMaximum(i)
        xiIndex(i) = xiIndex(i) - 1;
    end
    
    regionIndex = regionIndex + xiIndex(i)*len;
    len = len * (interpolant.xiLength(i)-1);
end

% translate coordinates to logical indices in hipercube
for i = 1:interpolant.inputDimensions
    xiHipercube(i) = (xi(i) - interpolant.xiMinimum(i))/interpolant.xiStep(i) - xiIndex(i);
end

% interpolate inside the region
[xo, membershipVector] = flhiInterpolateRegion(interpolant.regionVector(regionIndex).pointVector, interpolant.kernelFunctions, xiHipercube);

end

