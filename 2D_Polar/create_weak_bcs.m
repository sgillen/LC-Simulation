function [wbcMatrix] = create_weak_bcs(size)
    aMatrix = zeros(size);
    wMatrix = zeros(size);
    
    %aMatrix(:,end) = pi/4;
    %wMatrix(:,end) = 0;
    
    wbcMatrix(:,:,1) = aMatrix;
    wbcMatrix(:,:,2) = wMatrix;

end