function wbcMatrix = create_spherical_wbcs(size)

    aMatrix = zeros(size);
    bMatrix = zeros(size);
    wMatrix = zeros(size);
    
    aMatrix(:,:,1) = pi/2; 
    bMatrix(:,:,1) = pi/2;
    
    wMatrix(:,:,1) = 1;
    
    wbcMatrix(:,:,:,1) = aMatrix;
    wbcMatrix(:,:,:,2) = bMatrix;
    wbcMatrix(:,:,:,3) = wMatrix;
    
   


end