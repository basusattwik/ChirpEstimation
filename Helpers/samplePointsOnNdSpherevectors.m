function points = samplePointsOnNdSpherevectors(numPoints, radius, center)
    % num_points: number of points to sample
    % radius: radius of the hypersphere
    % center_point: an N-dimensional column vector representing the center of the hypersphere

    % Determine the dimension from the length of the center_point
    N = length(center);
    
    % Ensure center_point is a column vector
    center = center(:);

    % Generate points from a standard normal distribution
    randPoints= randn(N, numPoints); % N-dimensional points from a normal distribution

    % Normalize each point to lie on the surface of the unit hypersphere
    normFactor    = sqrt(sum(randPoints.^2, 1)); % Compute the norm of each point
    unitMagPoints = randPoints ./ normFactor;   % Normalize points to lie on the unit hypersphere

    % Scale the points to the desired radius
    scaledPoints = radius * unitMagPoints;

    % Shift the points to the desired center point
    points = bsxfun(@plus, scaledPoints, center);
end
