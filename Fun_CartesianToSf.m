%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project name: 
%Programer   : 
%Finish date :
%Records     :             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w,l_y,closestPoint] = Fun_CartesianToSf(point, path)

    A = path(1:end-1, :);
    B = path(2:end, :);

    BA = B - A;
    AP = point - A;
    dotProduct = dot(AP, BA, 2);
    normBA = vecnorm(BA, 2, 2).^2;
    t = dotProduct ./ normBA;
    t(t < 0) = 0;
    t(t > 1) = 1;
    projections = A + t .* BA;

    distances = vecnorm(projections - point, 2, 2);

    [minDistance, index] = min(distances);
    closestPoint = projections(index, :);

   if index < size(path, 1)
        directionVector = B(index, :) - A(index, :);
        normalVector = [-directionVector(2), directionVector(1)];  
        relativePosition = dot(point - closestPoint, normalVector);


         if relativePosition > 0
            l_y = -minDistance;
        elseif relativePosition < 0
            l_y = minDistance;
        else
            l_y = 0;
        end
   else
       l_y = 0;
   end
    
    w = (index);
end
