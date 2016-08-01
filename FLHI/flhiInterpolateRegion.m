function [xo, membershipVector] = flhiInterpolateRegion(pointVector, kernelFunctions, xi)
%Interpolates a position 'xi' in an enclosed space by points in 'pointVector'
%applying the associated kernel according to 'kernelFunctions'.
%   Something, something... dark side.
%
% pointVector.xi vector of input coordinates
% pointVector.xo vector of output coordinates
% 

inputDimensions = length(pointVector(1).xi);
outputDimensions = length(pointVector(1).xo);

xo = zeros(1, outputDimensions);

% pointVectorSize is 2^N, where N is the number of dimensions (inputDimensions)
pointVectorSize = length(pointVector);

membershipVector = zeros(1, pointVectorSize);

for xoIndex = 1:outputDimensions
    xoVector = zeros(pointVectorSize, 1);
    
    % complexity = (2^N)*N
    for pointIndex = 1:pointVectorSize
        point = pointVector(pointIndex);

        membership = zeros(1, inputDimensions);

        for xiIndex = 1:inputDimensions
            membership(xiIndex) = membershipEvaluate(point, kernelFunctions, xi, xiIndex, xoIndex);
        end

        membershipVector(pointIndex) = prod(membership);
        xoVector(pointIndex) = pointVector(pointIndex).xo(xoIndex);
    end
    
    % it isn't necessary to normalize membership functions since it's
    % assumed they have unitary sum from the rule "f + (1-f)"
    % in fact, this division might cause numeric errors
    %membershipVector = membershipVector./sum(membershipVector);

    xo(xoIndex) = membershipVector*xoVector;
end

end

function membership = kernelNearestNeighbor(xi)
    if(xi < 0.5)
        membership = 1;
    else
        membership = 0;
    end
end

function membership = kernelLinear(x)
    membership = 1 - x;
end

function membership = kernelCubic(x, f0, f1, df0, df1)
    % TODO: used to avoid numerical issues - is this really necessary?
    if(abs(f0-f1) < 1e-8)
        membership = 1;
    else
        membership = ((df0 + df1 + 2*f0 - 2*f1)*x^3 + (3*f1 - df1 - 3*f0 - 2*df0)*x^2 + df0*x + f0 - f1)/(f0 - f1);
    end
end

function membership = kernelLanczos(x)
    if x == 0
        membership = 1;
    elseif x > 0 && x < 1
        membership = 1*sin(pi*x)*sin(pi*x/1)/(pi^2*x^2);
    else
        membership = 0;
    end
end

function membership = kernelSpline(x)
    % from Keys, "Cubic Convolution Interpolation for Digital Image Processing," IEEE Transactions on Acoustics, Speech, and Signal Processing, Vol. ASSP-29, No. 6, December 1981, p. 1155. 
    membership = ((1.5 * x - 2.5) .* x) .* x + 1.0;
end

function membership = membershipEvaluate(point, kernelFunctions, xi, i, o)
    membership = 0;
    pxi = point.xi(i);
    kernelFunction = kernelFunctions(i);
    x = xi(i);
    
    switch kernelFunction
        case 0
            membership = kernelNearestNeighbor(x);
        case 1
            membership = kernelLinear(x);
        case 2
            membership = kernelCubic(x, point.cubicData(pxi+1).f0(i,o),  point.cubicData(pxi+1).f1(i,o),  point.cubicData(pxi+1).df0(i,o),  point.cubicData(pxi+1).df1(i,o));
        case 3
            membership = kernelLanczos(x);
        case 4
            membership = kernelSpline(x);
        otherwise
            error('Unknown membership function!');
    end
    
    % if the point is situated on the inverse side of the dimension, invert
    % the logic
    if(pxi == 1)
        membership = 1 - membership;
    end
    
    if(pxi ~= 0 && pxi ~= 1)
        error('Unknown base point position');
    end
end
