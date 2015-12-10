function [oMatrix] = enforce_2D_polar_SBCs(iMatrix)
%this is a pretty hacked together way of handling BC's but whatever. Just use comments to chose the BC's you want
%this returns oMatrix, which is the nMatrix with the strong boundaries
%enforced. Also returns seMatrix, which is a matrix containing the surface
%energy at specified boundary points. 

oMatrix = iMatrix;


%% Strong Boundaries
%end points
oMatrix(1,:) = -pi/4;
oMatrix(end,:) = pi/2;
%oMatrix(:,1) = pi/2;
%oMatrix(:,end)= pi/4;

end