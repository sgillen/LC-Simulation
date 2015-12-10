function oMatrix = enforce_2D_BCs(iMatrix)
%% this is a pretty hacked together way of handling BC's but whatever. Just use comments BC's to choose the BS's you want


%TODO, speed this up slightly by not copying iMatrix? can you even return
%values in matlab ? 
oMatrix = iMatrix ;

oMatrix(:,1,2) = 1;
oMatrix(:,1,1) = .5;
    
oMatrix(:,end,2) = -1;
oMatrix(:,end,1) = 0;

oMatrix(1,:,1) = 0;
oMatrix(1,:,2) = 1;

oMatrix(end,:,1)= 0;
oMatrix(end,:,2) = -1;

% nMatrix(:,1,2) = 0;
% nMatrix(:,end,2) = 0;
% nMatrix(1,1,1) = 1; nMatrix(1,1,2) = 1;
% nMatrix(end,1,1) = 1; nMatrix(end,1,2) = -1;
% nMatrix(1,end,1) = 1; nMatrix(1,end,2) = -1;
% nMatrix(end,end,1) = 1; nMatrix(end,end,2) = 1;


end