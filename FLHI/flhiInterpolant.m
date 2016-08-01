function [interpolant] = flhiInterpolant(pointVector, kernelFunctions)
%Receives a list of points in regular coordinates and returns a fuzzy
%logic hipercube interpolant
%   Something, something... dark side.

% TODO - validate input
%   - Verify if points are a regular (step is the same for all dimensions)
%   or semi-regular (different steps for dimensions), or irregular grid.
%   - Verify if points are enough to create regions (at least 2 points per
%   dimension)
%   - Verify if xiSteps are above zero
%   - Verify if xiSteps have equal rows (for regular and irregular grid)
%   - Verify if all points have equal xi and xo dimension sizes
%   - Verify if .xiLength is above one
%   - Verify if kernelFunctions length is one or inputDimensions
interpolant.pointVector = sortPoints(pointVector);

interpolant.pointCount = length(interpolant.pointVector);
interpolant.inputDimensions = length(interpolant.pointVector(1).xi);
interpolant.outputDimensions = length(interpolant.pointVector(1).xo);
interpolant.xiMinimum = interpolant.pointVector(1).xi;
interpolant.xiMaximum = interpolant.pointVector(end).xi;
interpolant.xiStep = countDimensionStep(interpolant.pointVector, interpolant.inputDimensions, interpolant.xiMinimum);
interpolant.xiLength = round((interpolant.xiMaximum - interpolant.xiMinimum)./interpolant.xiStep + 1);

% verify if xiLengths and pointCount are equal
if(interpolant.pointCount ~= prod(interpolant.xiLength))
    error('FLHI ERROR: Total number of input points does not match the expected point count');
end

% VERIFY
% if all(interpolant.xiStep == interpolant.xiStep(1)) REGULAR
% else, SEMIREGULAR
% TODO: IRREGULAR???

interpolant.regionCount = prod(interpolant.xiLength - 1);
%interpolant.regionVector = zeros(1, interpolant.regionCount);
interpolant.regionPointCount = 2^interpolant.inputDimensions;

if(length(kernelFunctions) == 1)
    kernelFunctions = kernelFunctions(1)*ones(1, interpolant.inputDimensions);
end

interpolant.kernelFunctions = kernelFunctions;

% calculate derivatives for all points (by central differences)
for i = 1:interpolant.pointCount
    interpolant.pointVector(i).dxo = zeros(interpolant.inputDimensions, interpolant.outputDimensions);
    
    for j = 1:interpolant.inputDimensions
        leftPointCoordinate = interpolant.pointVector(i).xi;
        rightPointCoordinate = interpolant.pointVector(i).xi;
        
        leftPointCoordinate(j) = leftPointCoordinate(j) - interpolant.xiStep(j);
        rightPointCoordinate(j) = rightPointCoordinate(j) + interpolant.xiStep(j);
        
        if(leftPointCoordinate(j) < interpolant.xiMinimum(j))
            leftPointCoordinate(j) = interpolant.xiMinimum(j);
        end
            
        if(rightPointCoordinate(j) > interpolant.xiMaximum(j))
            rightPointCoordinate(j) = interpolant.xiMaximum(j);
        end
        
        leftPointIndex = getIndexFromCoordinate(leftPointCoordinate, interpolant.inputDimensions, interpolant.xiMinimum, interpolant.xiStep, interpolant.xiLength);
        rightPointIndex = getIndexFromCoordinate(rightPointCoordinate, interpolant.inputDimensions, interpolant.xiMinimum, interpolant.xiStep, interpolant.xiLength);
        
        leftPoint = interpolant.pointVector(leftPointIndex);
        rightPoint = interpolant.pointVector(rightPointIndex);
        
        for k = 1:interpolant.outputDimensions
            %interpolant.pointVector(i).dxo(j,k) = (rightPoint.xo(k) - leftPoint.xo(k))/(2*interpolant.xiStep(j));
            interpolant.pointVector(i).dxo(j,k) = (rightPoint.xo(k) - leftPoint.xo(k))/2; % interpolation occurs at an unitary space which isn't dependent on step size
        end
    end
end

% add cubic data for all points
for i = 1:interpolant.pointCount
    interpolant.pointVector(i).cubicData(1).f0 = zeros(interpolant.inputDimensions, interpolant.outputDimensions);
    interpolant.pointVector(i).cubicData(1).f1 = zeros(interpolant.inputDimensions, interpolant.outputDimensions);
    interpolant.pointVector(i).cubicData(1).df0 = zeros(interpolant.inputDimensions, interpolant.outputDimensions);
    interpolant.pointVector(i).cubicData(1).df1 = zeros(interpolant.inputDimensions, interpolant.outputDimensions);
    interpolant.pointVector(i).cubicData(2).f0 = zeros(interpolant.inputDimensions, interpolant.outputDimensions);
    interpolant.pointVector(i).cubicData(2).f1 = zeros(interpolant.inputDimensions, interpolant.outputDimensions);
    interpolant.pointVector(i).cubicData(2).df0 = zeros(interpolant.inputDimensions, interpolant.outputDimensions);
    interpolant.pointVector(i).cubicData(2).df1 = zeros(interpolant.inputDimensions, interpolant.outputDimensions);

    for j = 1:interpolant.inputDimensions
        leftPointCoordinate = interpolant.pointVector(i).xi;
        rightPointCoordinate = interpolant.pointVector(i).xi;
        
        leftPointCoordinate(j) = leftPointCoordinate(j) - interpolant.xiStep(j);
        rightPointCoordinate(j) = rightPointCoordinate(j) + interpolant.xiStep(j);
        
        if(leftPointCoordinate(j) < interpolant.xiMinimum(j))
            leftPointCoordinate(j) = interpolant.xiMinimum(j);
        end
            
        if(rightPointCoordinate(j) > interpolant.xiMaximum(j))
            rightPointCoordinate(j) = interpolant.xiMaximum(j);
        end
        
        leftPointIndex = getIndexFromCoordinate(leftPointCoordinate, interpolant.inputDimensions, interpolant.xiMinimum, interpolant.xiStep, interpolant.xiLength);
        rightPointIndex = getIndexFromCoordinate(rightPointCoordinate, interpolant.inputDimensions, interpolant.xiMinimum, interpolant.xiStep, interpolant.xiLength);
        
        leftPoint = interpolant.pointVector(leftPointIndex);
        rightPoint = interpolant.pointVector(rightPointIndex);
        
        for k = 1:interpolant.outputDimensions
            interpolant.pointVector(i).cubicData(1).f0(j,k) = interpolant.pointVector(i).xo(k);
            interpolant.pointVector(i).cubicData(1).f1(j,k) = rightPoint.xo(k);
            interpolant.pointVector(i).cubicData(1).df0(j,k) = interpolant.pointVector(i).dxo(j,k);
            interpolant.pointVector(i).cubicData(1).df1(j,k) = rightPoint.dxo(j,k);
            interpolant.pointVector(i).cubicData(2).f0(j,k) = leftPoint.xo(k);
            interpolant.pointVector(i).cubicData(2).f1(j,k) = interpolant.pointVector(i).xo(k);
            interpolant.pointVector(i).cubicData(2).df0(j,k) = leftPoint.dxo(j,k);
            interpolant.pointVector(i).cubicData(2).df1(j,k) = interpolant.pointVector(i).dxo(j,k);
        end
    end
end


% TODO: this loop should be regionCount based, not pointCount
% Create regions
regionIndex = 1;
for i = 1:interpolant.pointCount
    basePoint = interpolant.pointVector(i);
    
    % if the base point is at the edges, this point can't form a region
    if(any(basePoint.xi == interpolant.xiMaximum))
        continue;
    end
    
    regionPointVector = interpolant.pointVector(i:i+interpolant.regionPointCount-1); % dirty hack to preallocate the structure vector
    
    pointShiftIndex = zeros(1, interpolant.inputDimensions);
    
    xoMatrix = zeros(interpolant.regionPointCount, interpolant.outputDimensions);
    
    for j = 1:interpolant.regionPointCount
        % Convert the point index and its shift to flat coordinates
        nextPointIndex = i; % starting from basePoint
        len = 1;
        for k = 1:interpolant.inputDimensions
            nextPointIndex = nextPointIndex + pointShiftIndex(k)*len;
            len = len * interpolant.xiLength(k);
        end
        
        % Add the shifted point to the region
        regionPointVector(j) = interpolant.pointVector(nextPointIndex);
        
        % Maintain shifted point in hipercube coordinates (0 and 1)
        regionPointVector(j).xi = round((regionPointVector(j).xi - basePoint.xi)./interpolant.xiStep);
        
        xoMatrix(j, :) = regionPointVector(j).xo;
        
        % Shift the index
        carry = 1;
        for k = 1:interpolant.inputDimensions
            pointShiftIndex(k) = pointShiftIndex(k) + carry;
            
            if(pointShiftIndex(k) > 1)
                pointShiftIndex(k) = 0;
                carry = 1;
            else
                break;
            end
        end
    end
    
    interpolant.regionVector(regionIndex).pointVector = regionPointVector;
    interpolant.regionVector(regionIndex).xoMaximum = max(xoMatrix);
    interpolant.regionVector(regionIndex).xoMinimum = min(xoMatrix);
    interpolant.regionVector(regionIndex).basePoint = basePoint;
    
    regionIndex = regionIndex + 1;
end

end

function [flatIndex] = getIndexFromCoordinate(coordinates, dimensionCount, minimum, step, length)

index = zeros(dimensionCount,1);

% translate coordinates to logical indices for table look up
for i = 1:dimensionCount
    % TODO: BOUNDS CHECK!
    index(i) = round((coordinates(i) - minimum(i))/step(i));
end

flatIndex = 1; % (start at 1, because Matlab...)
len = 1;
for i = 1:dimensionCount
    flatIndex = flatIndex + index(i)*len;
    len = len * length(i);
end


end

function [xiStep] = countDimensionStep(pointVector, inputDimensions, xiMinimum)

xiStep = zeros(1, inputDimensions);

pointCount = length(pointVector);

for i = 1:inputDimensions
    for j = 1:pointCount
        point = pointVector(j);
        if(point.xi(i) > xiMinimum(i))
            xiStep(i) = point.xi(i) - xiMinimum(i);
            break;
        end
    end
end

end