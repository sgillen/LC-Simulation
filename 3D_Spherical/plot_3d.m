function plot_3d(nMatrix)
    
    [X, Y, Z] = ndgrid(1:size(nMatrix,1),1:size(nMatrix,2),1:size(nMatrix,3));

    vt = squeeze(nMatrix(:,:,:,1));
    vp = squeeze(nMatrix(:,:,:,2));
    
    
    %in our case phi (vp) is the azithumal angle and theta (vp) our
    %elevation, which is actually the reverse of what matlab has
    figure
    [vx,vy,vz] = sph2cart(vp,vt,1);
    q3 = quiver3(X-vx./2,Y-vy./2,Z-vz./2,vx,vy,vz)
    
end

