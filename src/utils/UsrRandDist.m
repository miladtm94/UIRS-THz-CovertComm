%% Unifrmly distributed Random Points 
function [x, y]= UsrRandDist(xc,yc,rmin,rmax,numPoints)
 % Set your seed here if desired
 n = numPoints;
 r1 = rmin; r2 = rmax;
 r = sqrt(r1^2+(r2^2-r1^2)*rand(1,n)); % Using square root here ensures distribution uniformity by area
 for i=1:n
    t(i) = (2*pi/n)*i;
 end   
 x = xc + r.*cos(t); % Center of the circle in the x direction.
 y = yc + r.*sin(t); % Center of the circle in the y direction.
end