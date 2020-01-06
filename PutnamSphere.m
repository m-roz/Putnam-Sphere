% 1992 Putnam Exam Problem A-6:
% Four points are chosen at random on the surface of a sphere. What is the
% probability that the center of sphere lies inside the terahedron whose
% verticies are at the four points? (It is understood that each point is 
% independently chosen relative to a uniform distribution on the sphere)


% Using Monte Carlo simulation, this script computes the approximate
% probablity that a tetrahedron created from 4 random points on the surface
% of a sphere will encompass the origin of the sphere.

clear; clc;

% Number of trials
numTrials = 100000;

% For each trial, generate 4 random points on the surface of the sphere in
% spherical coordinates
r = 1;
phi = 2*pi*rand(4,numTrials);
theta = 2*pi*rand(4,numTrials);

% Convert to rectangular coordinates
x = r*sin(phi).*cos(theta);
y = r*sin(phi).*sin(theta);
z = r*cos(phi);

% Boolean vector for whether tetrahedron encompasses origin for each trial
origin = zeros(1,numTrials);

% Perform test for each trial
for i = 1:numTrials
    % Check if all points lie on any one side of the sphere. If they do, 
    % the origin must be outside the tetrahedron
    if all(x(:,i) > 0) | all(x(:,i) < 0) | all(y(:,i) > 0) | all(y(:,i) < 0) | all(z(:,i) > 0) | all(z(:,i) < 0)
        origin(i) = 0;
        
    else
    % Form a plane containing 3 of the points and find the perpendicular
    % distances between the plane and the origin and the plane and the 
    % 4th point. If the position vectors are both on the same side of the 
    % plane, the distances are compared
    
        % Form a vector VEC1 going from point 1 to point 2 that lies in the
        % plane
        VEC1 = [x(1,i) - x(2,i), y(1,i) - y(2,i), z(1,i) - z(2,i)];
        
        % Form a vector VEC2 going from point 1 to point 3 that lies in the
        % plane
        VEC2 = [x(1,i) - x(3,i), y(1,i) - y(3,i), z(1,i) - z(3,i)];
        
        % Take the cross product of the two vectors, producing a normal
        % vector N that establishes the orientation of the plane
        N = cross(VEC1,VEC2);
        
        % Find the perpendicular distance D_OriginToPlane from the 
        % spherical origin to the plane. This is done by taking a vector 
        % going from any point on the plane to the origin and then finding 
        % its projection onto vector N
        D_OriginToPlane = dot(N,[x(1,i), y(1,i), z(1,i)])/sqrt(dot(N,N));
        
        % Similarly, find the perpendicular distance D_Point4ToPlane from 
        % point 4 to the plane
        D_Point4ToPlane = dot(N,[x(4,i)-x(1,i), y(4,i)-y(1,i), z(4,i)-z(1,i)])/sqrt(dot(N,N));
        
        % Check the signs of the distances to see if both the origin
        % and point 4 are on the same side of the plane
        if sign(D_OriginToPlane) == sign(D_Point4ToPlane)
            % If the length from the plane to point 4 is greater than
            % the length from the plane to the origin, then the origin
            % must be contained within the tetrahedron
            origin(i) = abs(D_Point4ToPlane) >= abs(D_OriginToPlane);
        end
        
    end
end

% Calculate the approximate probability as a fraction of the number of 
% successful trials over the total number of trials
probability = sum(origin)/numTrials

% Actual probability is 0.1250 