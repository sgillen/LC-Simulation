function [wbcMatrix] = create_3d_wbcs(size)
    wbcMatrix = zeros(size(1),size(1),size(3),4);
    %We'll call wbcMatrix(:,:,:,1-3) the components of our weak vector and
    %wbcMatrix(:,:,:,4) the scaling factor
    
    %you'll need to normalize on your own (or don't, if you're into that)
    
    wbcMatrix(:,:,end,1) = sqrt(2);
    wbcMatrix(:,:,end,2) = sqrt(2);
    
    wbcMatrix(:,:,end,4) = sqrt(2);
    


end