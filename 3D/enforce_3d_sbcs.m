function oMatrix = enforce_3d_sbcs(iMatrix)
%% this is a pretty hacked together way of handling BC's but whatever. Just use comments to chose the BC's you want
%this returns oMatrix, which is the nMatrix with the strong boundaries
%enforced. Also returns seMatrix, which is a matrix containing the surface
%energy at specified boundary points. 
oMatrix = iMatrix;

oMatrix(:,:,1,1) = 0;
oMatrix(:,:,1,2) = 1;
oMatrix(:,:,1,3) = 0;

oMatrix(:,:,end,1) = .5;
oMatrix(:,:,end,2) = .5;
oMatrix(:,:,end,3) = 0;


end