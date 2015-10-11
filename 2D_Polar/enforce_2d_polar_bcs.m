function oMatrix = enforce_2D_polar_BCs(iMatrix)
%% this is a pretty hacked together way of handling BC's but whatever. Just use comments BC's to choose the BS's you want


%TODO, speed this up slightly by not copying iMatrix? can you even return

oMatrix = iMatrix;

oMatrix(:,1) = 0;

oMatrix(:,end) = pi/4;


end