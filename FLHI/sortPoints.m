function [pointVector] = sortPoints(pointVector)
%Sorts unordered points of 'pointVector', expecting a regular or
%semi-regular grid.
%   The return vector 'pointVector' contains points in a flat coordinate
%   system.

pointVectorSize = length(pointVector);

% Bubble sort (TODO: change this to a more efficient algorithm)
for i = 1:(pointVectorSize-1)
    for j = i+1:pointVectorSize
        if isPointGreater(pointVector(i), pointVector(j))
            % swap
            aux = pointVector(i);
            pointVector(i) = pointVector(j);
            pointVector(j) = aux;
        end
    end
end

end

function [result] = isPointGreater(p1, p2)
%Returns 1 if 'p1' is greater than 'p2', 0 otherwise.

result = 0;

%for i=1:length(p1.xi) % left to right order
for i=length(p1.xi):-1:1 % right to left order
    if p1.xi(i) > p2.xi(i)
        result = 1;
        return;
    elseif p1.xi(i) < p2.xi(i)
        result = 0;
        return;
    end
end

end
