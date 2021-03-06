%% initialize a grid and other parameters
%how long (in seconds) should each timestep be?
timestep = 0.0000001; 
%number of time steps to take
numsteps = 8000; 
%number of time frames to save
frames = 20;  

%how big is the grid (there is a director at every point)
grid = [5 5 5];
%how big is the total grid (meters?)
cellsize = [2e-6, 2e-6, 2e-6];
%spacing between directors
dx = cellsize(1)/grid(1);
dy = cellsize(2)/grid(2);
dz = cellsize(3)/grid(3);
%Creates a numbered grid thing
[X, Y, Z] = ndgrid(1:grid(1),1:grid(2),1:grid(3));

%bend, twist, and splay constants
K11 = 12e-12;
K22 = 5e-12;
K33 = 12e-12;

% rotational viscosity of 0.8 Poise
gamma = 0.08; 


%tolerance
tol = .001; 

convergent = false(1);
scaling_factor = .000001;



%% begin with some initial configuration

% each point in the matrix has an x y and z component, so nMatrix(x,y,z,1)
% is the x component of the point at (x,y,z) etc.

% We're about to hardcode so much stuff
 nMatrix = randn(grid(1)+2,grid(2)+2,grid(3),2);
 
 %Give an initial guess. 
 nMatrix(:,:,:,1) = 0;
 nMatrix(:,:,:,2) = 0;
 
 %nMatrix(:,:,1,:) = randn(grid(1)+2,grid(2)+2,2);


%% Create the boundary shape 


%should create a hemisphere shaped logical array, appears to work
%hemisphere = ((X-(grid(1)+1)/2).^2./((grid(1)-1)^2/4)+(Y-(grid(2)+1)/2).^2./((grid(2)-1)^2/4)+(Z./grid(3)).^2) < 1;

%should be an octant of a sphere, this is a good place to start if things
%look like they're not working though
octant = sqrt(((X-1)/(grid(1)-1)).^2 + ((Y-1)/(grid(2)-1)).^2 + ((Z-1)/(grid(3)-1)).^2) < 1;

%initialize edge, which is the outermost layer of the hemisphere
bclocation(:,:,:,1) = (~octant);
bclocation(:,:,:,2) = (~octant);

%a hacked together but usable way to get the outermost of (I think) any
%shape, should probaly throw it into a function

%{
for ii = 1:(size(octant,1)-1)
   for jj = 1:(size(octant,2)-1)
       for k = 1:size(octant,3)
           
           if(octant(ii,jj,k) == 1) %only pick out points which were previously inside the octant                                                 
                if(octant(ii,jj+1,k) == 0)                                                                                                                
                     edge(ii,jj,k,:) = 1;
                end
                if (octant(ii+1,jj,k) == 0)
                     edge(ii,jj,k,:) = 1;
                end  
           end 
           
       end
   end
end
%}

%{

%you can use this for the hemisphere if you want
%{
for ii = 2:(size(octant,1)-1)
   for jj = 2:(size(octant,2)-1)
       for k = 1:size(octant,3)
           
           if(octant(ii,jj,k) == 1) %only pick out points which were previously inside the octant                                                 
                if(octant(ii,jj+1,k) + octant(ii,jj-1,k) == 1)                                                                                                                
                     edge(ii,jj,k,:) = 1;
                elseif (octant(ii+1,jj,k) + octant(ii-1,jj,k) == 1)
                     edge(ii,jj,k,:) = 1;
                end
                
           end 
           
       end
   end
end
%}



%This is Davids way of getting the edge, seems to produce the same
%results, I'll leave it here for now, but use mine since I understand it
%better

%Davids way also doesn't quite work for the octant
%{
Nz = grid(3);
Nx = grid(1);
Ny = grid(2);

for zz=(1:Nz-1)

    B(zz)=bwboundaries(octant(:,:,zz),'noholes');

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
%}


theta = zeros(grid(1),grid(2), grid(3));
phi   = zeros(grid(1),grid(2), grid(3));
flat  = ~octant(:,:,1);

%we need to use atan2 here, which is the 4 quadrent arctan. normally we use
%theta = atan(y/x) with atan2 it's theta = atan2(y,x). this allows matlab
%to vary the angle between 0 and 2pi rather than -pi/2 to pi/2


%suspect phi is wrong
%atan2((Y-(grid(2)+1)/2),(X-(grid(1)+1)/2));
phi         = phi+bclocation(:,:,:,1).*atan2(Y-1,X-1);

%acos((Z-(grid(3)+1)/2)./sqrt((X-(grid(1)+1)/2).^2+(Y-(grid(2)+1)/2).^2+(Z-(grid(3)+1)/2).^2));
theta       = theta+bclocation(:,:,:,2).*atan2(Z-1, sqrt((X-1).^2 + (Y-1).^2));
theta(flat) = 0;

bcangles(:,:,:,1) = theta;
bcangles(:,:,:,2) = phi;

%had some problems with NaN, we can just set them to zero though
%edge(isnan(edge)) = 0;    


%array of normal vectors, of course it only really makes sense to define a
%normal vector on the surface, but since we're only using those points
%anyway that won't really matter
%norm = zeros(grid(1), grid(2), grid(3), 3);
%norm(:,:,:,1) = sin(theta).*cos(phi).*hemisphere;   %.*((xscaled>0)-0.5)*2;   %David had this line here at the end,I'm not too sure we'll need it, let's see
%norm(:,:,:,2) = sin(theta).*sin(phi).*hemisphere;  %.*((xscaled>0)-0.5)*2;    
%norm(:,:,:,3) = cos(theta).*hemisphere;


%% set some boundary conditions. these should be enforced at each step.

nMatrix = enforce_3d_spherical_sbcs(nMatrix,bcangles,bclocation);


s(1) = size(nMatrix,1);
s(2) = size(nMatrix,2);
s(3) = size(nMatrix,3);
wbcMatrix = create_spherical_wbcs(s);

nMatrix = mod(nMatrix,2*pi); 
%% !!! if you want to call a mod make sure you do it before the unwraping, otherwise you defeat the purpose. 
% unwrap our angle, this should be done at every step
nMatrix(:,:,:,1) = unwrap(nMatrix(:,:,:,1),pi,1);
nMatrix(:,:,:,1) = unwrap(nMatrix(:,:,:,1),pi,2);
nMatrix(:,:,:,1) = unwrap(nMatrix(:,:,:,1),pi,3);


nMatrix(:,:,:,2) = unwrap(nMatrix(:,:,:,2),pi,1);
nMatrix(:,:,:,2) = unwrap(nMatrix(:,:,:,2),pi,2);
nMatrix(:,:,:,2) = unwrap(nMatrix(:,:,:,2),pi,3);

%%


% initialize arrays that we can use to watch the system evolve
vs = repmat(zeros(size(nMatrix)),1,1,1,1,frames);
vs(:,:,:,:,1) = nMatrix;

% take a picture whenever we hit one of the following time periods
maxstep = zeros(1,frames);
snapshot = round(linspace(1,numsteps,frames))';
snapnum = 1;

avgEnergy = zeros(1,frames);
avgEnergy(snapnum) = mean2(lc_energy_spherical(nMatrix, K11, K22, K33, dx, dy, dz));

snapnum = 2;

%% handle the calculations
for ii = 2:numsteps
    % calulate step sizes
    %grab the EL terms
    [ELt, ELp]= EL_terms_spherical(nMatrix,K11,K22,K33,dx,dy,dz,wbcMatrix);
    nMatrix(:,:,:,1) = nMatrix(:,:,:,1) - timestep/gamma*ELt;
    nMatrix(:,:,:,2) = nMatrix(:,:,:,2) - timestep/gamma*ELp;

    
    %% enforce BC's
    %should maybe make enforcing BCs a function. 
    %nMatrix(:,:,1,1) = Y(:,:,1) - (grid(2)+1)/2;
    %nMatrix(:,:,1,2) = -X(:,:,1) + (grid(1)+1)/2;
    %nMatrix(:,:,1,3) = 0;
    
    %nMatrix(:,:,:,1) = nMatrix(:,:,:,1).*(hemisphere);
    %nMatrix(:,:,:,2) = nMatrix(:,:,:,2).*(hemisphere);
    %nMatrix(:,:,:,3) = nMatrix(:,:,:,3).*(hemisphere);
    
    nMatrix = enforce_3d_spherical_sbcs(nMatrix,bcangles,bclocation);
    
    %for sure find a way to vectorize this. 
    %nMatrix = setstrongboundary(nMatrix, edge, norm);
   
    % renormalize and replace NaN with 0
    nMatrix(isnan(nMatrix)) = 0;  
    nMatrix = mod(nMatrix,2*pi);
    
    %unwrap our angles. This actually takes a LOT of time, just fyi.
    nMatrix(:,:,:,1) = unwrap(nMatrix(:,:,:,1),pi,1);
    nMatrix(:,:,:,1) = unwrap(nMatrix(:,:,:,1),pi,2);
    nMatrix(:,:,:,1) = unwrap(nMatrix(:,:,:,1),pi,3);


    nMatrix(:,:,:,2) = unwrap(nMatrix(:,:,:,2),pi,1);
    nMatrix(:,:,:,2) = unwrap(nMatrix(:,:,:,2),pi,2);
    nMatrix(:,:,:,2) = unwrap(nMatrix(:,:,:,2),pi,3);
   
   
    
    
    %save at the certain points we care about. 
    if any(ii==snapshot)      
       avgEnergy(snapnum) = mean2(lc_energy_spherical(nMatrix, K11, K22, K33, dx, dy, dz));
        fprintf('step %d/%d      avg = %d  timestep = %d\n', snapnum, frames, avgEnergy(snapnum), timestep);
        
        %try this when you're feeling cheeky
        %{
          if(~convergent)
               if(abs(avgEnergy(snapnum) - avgEnergy(snapnum - 1))) < tol
                   convergent = true(1);
                   %timestep = .1;
                   test_value = avgEnergy(snapnum);
                   fprintf('things look convergent at  a value of : %d, trying larger timestep\n', avgEnergy(snapnum))
               else
                   timestep = abs(scaling_factor*(avgEnergy(snapnum) - avgEnergy(snapnum - 1)));
               end
            
    elseif(convergent)
          if(avgEnergy(snapnum) - avgEnergy(snapnum - 1) < tol)
               if(abs(test_value - avgEnergy(snapnum)) < tol)
                   fprintf('looks like we are within the tolerance')
                   break % break out of the outer while loop  
                else
                   fprintf('possible convergent value found %d \n', avg(snapnum));
                   timestep = .1;
                   test_value = avgEnergy(snapnum); 
                end
          else
                 timestep = abs(scaling_factor*(avgEnergy(snapnum) - avgEnergy(snapnum - 1)));
          end
           
     end
    %}
        vs(:,:,:,:,snapnum) = nMatrix;
        snapnum = snapnum + 1;
    end
end

%% analysis
%vtAll = squeeze(vs(:,:,:,1,:));
%vpAll = squeeze(vs(:,:,:,2,:));

%vt = squeeze(vtAll(:,:,:,end));
%vp = squeeze(vpAll(:,:,:,end));

%vx = sin(vp).*sin(vt);
%vy = cos(vp).*sin(vt);
%vz = cos(vt);

%figure(1)
%clf
%slice = 1;
%q2 = quiver(X(:,:,slice)-vx(:,:,slice)./2,Y(:,:,slice)-vy(:,:,slice)./2,vx(:,:,slice),vy(:,:,slice));
%set(q2,'ShowArrowHead','off')

figure(2)
clf
%I wrote this one myself, make sure it's in matlabs path.
plot_3d(nMatrix);

%[vx,vy,vz] = sph2cart(vt,vp,1);
%questionable whether this is really working 
%q3 = quiver3(X-vx(3:end)./2,Y-vy(3:end)./2,Z-vz(3:end)./2,vx,vy,vz);
%set(q3,'ShowArrowHead','off')

figure(3)
clf
plot(avgEnergy(2:end));


%% let's throw this to mathematica for more analysis
save('/home/gillen/Documents/Computation/3D_Spherical/n3d.mat','vs');
fprintf('saved\n');