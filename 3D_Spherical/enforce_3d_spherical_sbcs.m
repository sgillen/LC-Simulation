function oMatrix = enforce_3d_spherical_sbcs(iMatrix,angles,location)
%% this is a pretty hacked together way of handling BC's but whatever. Just use comments BC's to choose the BC's you want
%TODO, speed this up slightly by not copying iMatrix? can you even return

oMatrix = iMatrix;


%Bottom strong BC, enforce this frst if you're going to and then let the
%edge take care of the rest. 
oMatrix(:,:,1,1) = 0;
oMatrix(:,:,1,2) = pi;


tmpMatrix = oMatrix(3:end,3:end,:,:);
tmpMatrix((location ~= 0)) = angles((location ~= 0));
oMatrix(3:end,3:end,:,:) = tmpMatrix;

%hardcoded
%theta 
oMatrix(2,:,:,1) = (pi - oMatrix(3,:,:,1));
oMatrix(1,:,:,1) = (pi - oMatrix(4,:,:,1));

%phi
oMatrix(2,:,:,2) = 2*pi - oMatrix(3,:,:,2);
oMatrix(1,:,:,2) = 2*pi - oMatrix(4,:,:,2);


%theta 
oMatrix(:,2,:,1) =  oMatrix(:,3,:,1);
oMatrix(:,1,:,1) =  oMatrix(:,4,:,1);

%phi
oMatrix(:,2,:,2) = -oMatrix(:,3,:,2);
oMatrix(:,1,:,2) = -oMatrix(:,4,:,2);



%oMatrix(:,:,1,1) = pi/4;
%oMatrix(:,:,1,2) = pi/4;

%oMatrix(:,:,end,1) = pi/4;
%oMatrix(:,:,end,2) = pi/4;


%oMatrix(:,:,end,1) = pi;
%oMatrix(:,:,end,2) = pi;

%oMatrix(:,1,:,1) = pi;
%oMatrix(:,1,:,2) = pi;

%oMatrix(:,end,:,1) = pi;
%oMatrix(:,end,:,2) = pi;

%oMatrix(1,:,:,1) = pi;
%oMatrix(1,:,:,2) = pi;

%oMatrix(end,:,:,1) = pi;
%oMatrix(end,:,:,2) = pi;



end