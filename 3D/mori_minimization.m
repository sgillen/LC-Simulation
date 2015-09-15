%%The majority of this code was written by David Somers, for now I'm
%merely making some modifications

%% initialize a grid and other parameters
numsteps = 2000; %steps in time to take
frames = 20;  %number of time frames to save
grid = [12 12 8]; %how big is the grid (there is a director at every point)
timestep = 0.0005; %how long (in seconds) should each timestep be?
gamma = 0.08; % rotational viscosity of 0.8 Poise
cellsize = [20e-6, 20e-6, 10e-6];
dx = cellsize(1)/grid(1);
dy = cellsize(2)/grid(2);
dz = cellsize(3)/grid(3);
K11 = 12e-12;
K22 = 5e-12;
K33 = 12e-12;

[X, Y, Z] = ndgrid(1:grid(1),1:grid(2),1:grid(3));

hemisphere = ((X-(grid(1)+1)/2).^2./((grid(1)-1)^2/4)+(Y-(grid(2)+1)/2).^2./((grid(2)-1)^2/4)+(Z./grid(3)).^2) < 1;

 
boundary = zeros(grid(1), grid(2), grid(3));
for ii = 2:(size(hemisphere,1)-1)
   for jj = 2:(size(hemisphere,2)-1)
       for k = 1:size(hemisphere,3)
           
           if(hemisphere(ii,jj,k) == 1) %only pick out points which were previously inside the hemisphere                                                 
                if(hemisphere(ii,jj+1,k) + hemisphere(ii,jj-1,k) == 1)                                                                                                                
                     boundary(ii,jj,k) = 1;
                elseif (hemisphere(ii+1,jj,k) + hemisphere(ii-1,jj,k) == 1)
                     boundary(ii,jj,k) = 1;
                end
                
           end 
           
       end
   end
end


theta = zeros(grid(1),grid(2), grid(3));
phi   = zeros(grid(1),grid(2), grid(3));

phi         = phi+hemisphere.*atan((Y-(Ny+1)/2)./(X-(Nx+1)/2));
theta       = theta+hemisphere.*acos((Z-(Nz+1)/2)./sqrt((X-(Nx+1)/2).^2+(Y-(Ny+1)/2).^2+(Z-(Nz+1)/2).^2));
theta(flat) = pi/2;
    

%array of normal vectors, of course it only really makes sense to define a
%normal vector on the surface, but since we're only using those points
%anyway that won't really matter
norm = zeros(grid(1), grid(2), grid(3), 3);

norm(:,:,:,1) = sin(theta).*cos(phi).*hemisphere;   %.*((xscaled>0)-0.5)*2;   %David had this line here at the end,I'm not too sure we'll need it, let's see
norm(:,:,:,2) = sin(theta).*sin(phi).*hemisphere;  %.*((xscaled>0)-0.5)*2;    
norm(:,:,:,3) = cos(theta).*hemisphere;


%{
This is Davids way of getting the boundary, seems to produce the same
results, I'll leave it here for now, but use mine since I understand it
better

Nz = grid(3);
Nx = grid(1);
Ny = grid(2);

for zz=(1:Nz-1)

    B(zz)=bwboundaries(hemisphere(:,:,zz),'noholes');

end

 

edge=zeros(Nx,Ny,Nz);

for zz=(1:Nz-1)

    for ii=1:length(B{zz})

        edgex=B{zz}(ii,1);

        edgey=B{zz}(ii,2);

        edge(edgex,edgey,zz)=1;

    end

end

%}
 
%boundary
%edge(:,:,end-1)=hemisphere(:,:,end-1)

%% begin with some initial configuration

% each point in the matrix has an x y and z component, so nMatrix(x,y,z,1)
% is the x component of the point at (x,y,z) etc.


% nMatrix = ones(grid(1),grid(2),grid(3),3);
% nMatrix(:,:,:,1) = 0;
% nMatrix(:,:,:,3) = 0;

% nMatrix(:,:,:,1) = repmat(Y(:,:,1) - (grid(2)+1)/2,1,1,grid(3));
% nMatrix(:,:,:,2) = repmat(-X(:,:,1) + (grid(1)+1)/2,1,1,grid(3));
% nMatrix(:,:,:,3) = 0;


 nMatrix = randn(grid(1),grid(2),grid(3),3);
 nMatrix(:,:,:,3) = 0;

%% set some boundary conditions. these should be enforced at each step.
%nMatrix(:,:,1,1) = 1; %Y(:,:,1) - (grid(2)+1)/2;
%nMatrix(:,:,1,2) = 0; %-X(:,:,1) + (grid(1)+1)/2;
%nMatrix(:,:,1,3) = 0;

%nMatrix(:,:,end,1) = 0;
%nMatrix(:,:,end,2) = 1;
%nMatrix(:,:,end,3) = 0;

%kind of like a very strong torque
nMatrix(:,:,1,1) = 0;
nMatrix(:,:,1,2) = 1;
nMatrix(:,:,1,3) = 0;

%hemisphere is zero outside the boundary so this effectively sets
%everything outside the droplet to zero 

nMatrix(:,:,:,1) = nMatrix(:,:,:,1).*double(hemisphere);
nMatrix(:,:,:,2) = nMatrix(:,:,:,2).*double(hemisphere);
nMatrix(:,:,:,3) = nMatrix(:,:,:,3).*double(hemisphere);


%nMatrix = setBoundary(nMatrix, boundary, grad);


%% do some other things
% normalize. this should also be done at each step.

%getting NaN when dividing by zero, the little logical indexing bit below
%is a workaround for that
nMatrix = nMatrix./repmat(sqrt(sum(nMatrix.^2,4)),[1 1 1 3]);
nMatrix(isnan(nMatrix)) = 0;  

% initialize arrays that we can use to watch the system evolve
vs = repmat(zeros(size(nMatrix)),1,1,1,1,frames);


vs(:,:,:,:,1) = nMatrix;

% take a picture whenever we hit one of the following time periods
maxstep = zeros(1,frames);
snapshot = round(linspace(1,numsteps,frames))';
snapnum = 2;



%% handle the calculations
for ii = 2:numsteps
    % calulate step sizes
    [ELx, ELy, ELz] = EL_terms(nMatrix,K11,K22,K33,dx,dy,dz);
    nMatrix(:,:,:,1) = nMatrix(:,:,:,1) - timestep/gamma*ELx;
    nMatrix(:,:,:,2) = nMatrix(:,:,:,2) - timestep/gamma*ELy;
    nMatrix(:,:,:,3) = nMatrix(:,:,:,3) - timestep/gamma*ELz;
    
    % enforce BC's
    %nMatrix(:,:,1,1) = Y(:,:,1) - (grid(2)+1)/2;
    %nMatrix(:,:,1,2) = -X(:,:,1) + (grid(1)+1)/2;
    %nMatrix(:,:,1,3) = 0;
    
    
    %nMatrix(:,:,end,1) = 0;
    %nMatrix(:,:,end,2) = 1;
    %nMatrix(:,:,end,3) = 0;
    
    nMatrix(:,:,1,1) = 0;   
    nMatrix(:,:,1,2) = 1;
    nMatrix(:,:,1,3) = 0;

    
    nMatrix(:,:,:,1) = nMatrix(:,:,:,1).*(hemisphere);
    nMatrix(:,:,:,2) = nMatrix(:,:,:,2).*(hemisphere);
    nMatrix(:,:,:,3) = nMatrix(:,:,:,3).*(hemisphere);
    
   nMatrix = setBoundary(nMatrix, boundary, norm);


   
    % renormalize and replace NaN with 0
    nMatrix = nMatrix./repmat(sqrt(sum(nMatrix.^2,4)),[1 1 1 3]);
    nMatrix(isnan(nMatrix)) = 0;  
    
    if any(ii==snapshot)
        fprintf('%d%%\n',round(100*ii/numsteps))
        maxstep(snapnum)=timestep/gamma*max([max(ELx(:)),max(ELy(:)),...
            max(ELz(:))]);
        fprintf('%d\n',maxstep(snapnum))
        vs(:,:,:,:,snapnum) = nMatrix;
        snapnum = snapnum + 1;
    end
end

%% analysis
vxAll = squeeze(vs(:,:,:,1,:));
vyAll = squeeze(vs(:,:,:,2,:));
vzAll = squeeze(vs(:,:,:,3,:));

vx = squeeze(vxAll(:,:,:,end));
vy = squeeze(vyAll(:,:,:,end));
vz = squeeze(vzAll(:,:,:,end));

figure(1)
clf
slice = 1;
q2 = quiver(X(:,:,slice)-vx(:,:,slice)./2,Y(:,:,slice)-vy(:,:,slice)./2,vx(:,:,slice),vy(:,:,slice));
set(q2,'ShowArrowHead','off')

figure(2)
clf
q3 = quiver3(X-vx./2,Y-vy./2,Z-vz./2,vx,vy,vz);
set(q3,'ShowArrowHead','off')

figure(3)
clf
plot(maxstep(2:end))


%% let's throw this to mathematica for more analysis
save('/home/gillen/Documents/Computation/3D/n.mat','vs');
fprintf('saved\n');
%{
% % this will find the jones transmission matrix for each point in (x,y)
% jones_matrices = LC_matrices(nMatrix);
% % for convenience, reshape the (x,y,2x2) matrix into (2x2,x,y). This will
% % save some squeezing later.
% jones_matrices = permute(jones_matrices,[3,4,1,2]);
% input_ray = [1;0];
% output_ray = zeros(2,1,grid(1),grid(2));
% for xx = 1:grid(1)
%     for yy = 1:grid(2)
%         output_ray(:,:,xx,yy) = jones_matrices(:,:,xx,yy)*input_ray;
%     end
% end
%
%
% [Xq, Yq] = meshgrid(linspace(-1,1),linspace(-1,1));
% % Iq = interp2(X, Y, cos(squeeze(output_ray(1,:,:,:))).^2, Xq, Yq);
%}