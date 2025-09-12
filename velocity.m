function [ velocity ] = velocity( currentPosition, goalPosition )
%VELOCITY Summary of this function goes here
%   Detailed explanation goes here
    paramMaxVelocity = 2;
    paramReachDistance = 4;
    paramSmoothFactor = 2;
    velocity = [0; 0; 0];

    errorPosition = goalPosition - currentPosition;

    errorNorm = sqrt(sum(errorPosition.^2));
    
    if errorNorm > 1.0
       
        scale = paramMaxVelocity/errorNorm;
        velocity= errorPosition.*scale;

        if errorNorm < paramReachDistance
            velocity = velocity *  ((errorNorm./paramReachDistance).^paramSmoothFactor);
        end
    end

end

