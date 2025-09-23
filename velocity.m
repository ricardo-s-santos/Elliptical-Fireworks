function [ velocity ] = velocity( currentPosition, goalPosition)
    paramMaxVelocity = 0.1;
    paramReachDistance = 0.4;
    paramSmoothFactor = 0.1;
    velocity = [0; 0; 0];

    errorPosition = goalPosition - currentPosition;

    errorNorm = sqrt(sum(errorPosition.^2));
    
    if errorNorm > 0.1     
        scale = paramMaxVelocity/errorNorm;
        velocity= errorPosition.*scale;
        if errorNorm < paramReachDistance
            velocity = velocity *  ((errorNorm./paramReachDistance).^paramSmoothFactor);
        end
    end

end

