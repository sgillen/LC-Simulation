function [ bclocation bcangles ] = create_boundary_shape()
% returns two matrices representing which locations have a strong boundary
% and the direction of that boundary

%as with everything else in this code base, just use comments to select
%which shape you want 



%octant
%{
octant = sqrt(((X-1)/(grid(1)-1)).^2 + ((Y-1)/(grid(2)-1)).^2 + ((Z-1)/(grid(3)-1)).^2) < 1;

bclocation = zeros(size(octant));
bclocation(:,:,:,1) = (~octant);
bclocation(:,:,:,2) = (~octant);

theta = zeros(grid(1),grid(2), grid(3));
phi   = zeros(grid(1),grid(2), grid(3));
flat  = ~octant(:,:,1);


%we need to use atan2 here, which is the 4 quadrent arctan. normally we use
%theta = atan(y/x) with atan2 it's theta = atan2(y,x). this allows matlab
%to vary the angle between 0 and 2pi rather than -pi/2 to pi/2

phi         = phi+bclocation(:,:,:,1).*atan2(Y-1,X-1);
theta       = theta+bclocation(:,:,:,2).*atan2(Z-1, sqrt((X-1).^2 + (Y-1).^2));
theta(flat) = 0;

bcangles(:,:,:,1) = theta;
bcangles(:,:,:,2) = phi;
%}

%hemispshere
hemisphere = ((X-(grid(1)+1)/2).^2./((grid(1)-1)^2/4)+(Y-(grid(2)+1)/2).^2./((grid(2)-1)^2/4)+(Z./grid(3)).^2) < 1;

bclocation = zeros(size(hemisphere));
bclocation(:,:,:,1) = (~hemisphere);
bclocation(:,:,:,2) = (~hemisphere);


%we need to use atan2 here, which is the 4 quadrent arctan. normally we use
%theta = atan(y/x) with atan2 it's theta = atan2(y,x). this allows matlab
%to vary the angle between 0 and 2pi rather than -pi/2 to pi/2

phi = phi+bclocation(:,:,:,1).*atan2((Y-(grid(2)+1)/2),(X-(grid(1)+1)/2));

theta = theta+bclocation(:,:,:,2).*acos((Z-(grid(3)+1)/2)./sqrt((X-(grid(1)+1)/2).^2+(Y-(grid(2)+1)/2).^2+(Z-(grid(3)+1)/2).^2));
theta(flat) = 0;

bcangles(:,:,:,1) = theta;
bcangles(:,:,:,2) = phi;

end

