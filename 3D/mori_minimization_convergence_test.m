%%The majority of this code was written by David Somers, for now I'm
%merely making some modifications

%% initialize a grid and other parameters
%how long (in seconds) should each timestep be?
timestep = 0.0025; 
%number of time steps to take
numsteps = 4000; 
%number of time frames to save
frames = 20;  

%how big is the grid (there is a director at every point)
grid = [12  12  8];
%how big is the total grid (meters?)
cellsize = [20e-6, 20e-6, 10e-6];
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
tol = .01; 

convergent = false(1);
scaling_factor = 1;

%% begin with some initial configuration

% each point in the matrix has an x y and z component, so nMatrix(x,y,z,1)
% is the x component of the point at (x,y,z) etc.

 nMatrix = randn(grid(1),grid(2),grid(3),3);
 %nMatrix(:,:,:,3) = 0;


%% Create the boundary shape 
%{
%should create a hemisphere shaped logical array, appears to work
hemisphere = ((X-(grid(1)+1)/2).^2./((grid(1)-1)^2/4)+(Y-(grid(2)+1)/2).^2./((grid(2)-1)^2/4)+(Z./grid(3)).^2) < 1;

%initialize boundary, which is the outermost layer of the hemisphere
boundary(:,:,:,1) = zeros(grid(1), grid(2), grid(3));
boundary(:,:,:,2) = zeros(grid(1), grid(2), grid(3));
boundary(:,:,:,3) = zeros(grid(1), grid(2), grid(3));

%a hacked together but usable way to get the outermost of (I think) any
%shape, should probaly throw it into a function
for ii = 2:(size(hemisphere,1)-1)
   for jj = 2:(size(hemisphere,2)-1)
       for k = 1:size(hemisphere,3)
           
           if(hemisphere(ii,jj,k) == 1) %only pick out points which were previously inside the hemisphere                                                 
                if(hemisphere(ii,jj+1,k) + hemisphere(ii,jj-1,k) == 1)                                                                                                                
                     boundary(ii,jj,k,:) = 1;
                elseif (hemisphere(ii+1,jj,k) + hemisphere(ii-1,jj,k) == 1)
                     boundary(ii,jj,k,:) = 1;
                end
                
           end 
           
       end
   end
end

theta = zeros(grid(1),grid(2), grid(3));
phi   = zeros(grid(1),grid(2), grid(3));
flat  = hemisphere(:,:,1);

%we need to use atan2 here, which is the 4 quadrent arctan. normally we use
%theta = atan(y/x) with atan2 it's theta = atan2(y,x). this allows matlab
%to vary the angle between 0 and 2pi rather than -pi/2 to pi/2
phi         = phi+hemisphere.*atan2((Y-(grid(2)+1)/2),(X-(grid(1)+1)/2));
theta       = theta+hemisphere.*acos((Z-(grid(3)+1)/2)./sqrt((X-(grid(1)+1)/2).^2+(Y-(grid(2)+1)/2).^2+(Z-(grid(3)+1)/2).^2));
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
 
%}
%% set some boundary conditions. these should be enforced at each step.

%kind of like a very strong torque
nMatrix(:,:,1,1) = 0;
nMatrix(:,:,1,2) = 1;
nMatrix(:,:,1,3) = 0;

nMatrix(:,:,end,1) = .5;
nMatrix(:,:,end,2) = .5;
nMatrix(:,:,end,3) = 0;


%hemisphere is zero outside the boundary so this effectively sets
%everything outside the droplet to zero 

%{
nMatrix(:,:,:,1) = nMatrix(:,:,:,1).*double(hemisphere);
nMatrix(:,:,:,2) = nMatrix(:,:,:,2).*double(hemisphere);
nMatrix(:,:,:,3) = nMatrix(:,:,:,3).*double(hemisphere);

nMatrix = setstrongboundary(nMatrix, boundary, norm);
%}




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
    %grab the EL terms
    [ELx, ELy, ELz] = EL_terms(nMatrix,K11,K22,K33,dx,dy,dz);
    nMatrix(:,:,:,1) = nMatrix(:,:,:,1) - timestep/gamma*ELx;
    nMatrix(:,:,:,2) = nMatrix(:,:,:,2) - timestep/gamma*ELy;
    nMatrix(:,:,:,3) = nMatrix(:,:,:,3) - timestep/gamma*ELz;
    
    %% enforce BC's
    %should maybe make enforcing BCs a function. 
    %nMatrix(:,:,1,1) = Y(:,:,1) - (grid(2)+1)/2;
    %nMatrix(:,:,1,2) = -X(:,:,1) + (grid(1)+1)/2;
    %nMatrix(:,:,1,3) = 0;
    
     nMatrix(:,:,1,1) = 0;
     nMatrix(:,:,1,2) = 1;
     nMatrix(:,:,1,3) = 0;

     nMatrix(:,:,end,1) = .5;    
     nMatrix(:,:,end,2) = .5;
     nMatrix(:,:,end,3) = 0;
    
    %nMatrix(:,:,:,1) = nMatrix(:,:,:,1).*(hemisphere);
    %nMatrix(:,:,:,2) = nMatrix(:,:,:,2).*(hemisphere);
    %nMatrix(:,:,:,3) = nMatrix(:,:,:,3).*(hemisphere);
    
    
   %for sure find a way to vectorize this. 
   %nMatrix = setstrongboundary(nMatrix, boundary, norm);
   
    % renormalize and replace NaN with 0
    nMatrix = nMatrix./repmat(sqrt(sum(nMatrix.^2,4)),[1 1 1 3]);
    nMatrix(isnan(nMatrix)) = 0;  
    
    %save at the certain points we care about. 
    if any(ii==snapshot)
        fprintf('%d%%\n',round(100*ii/numsteps))
        avg(snapnum) = mean2(lc_energy(nMatrix, K11, K22, K33, dx, dy, dz));
        fprintf('%d\n', avg(snapnum))
        
        if(~convergent)
            if(abs(avg(snapnum) - avg(snapnum - 1))) < tol
                convergent = true(1);
                timestep = .1;
                test_value = avg(snapnum);
            else
                timestep = 1*avg(snapnum);
            end
            
        elseif(convergent)
            fprintf('things are convergent\n')
             if(avg(snapnum) - avg(snapnum - 1) < tol)
                 if(abs(test_value - avg(snapnum)) < tol)
                    fprintf('hello! we should be breaking now\n')
                    break % break out of the 
                 else
                    test_value = avg(snapnum) 
                 end
             end
           
        end
    
        vs(:,:,:,:,snapnum) = nMatrix;
        snapnum = snapnum + 1;
    end
end


%% TODO : make analys work for n length arrays

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
plot(avg(2:end))


%% let's throw this to mathematica for more analysis
save('/home/gillen/Documents/Computation/3D/n.mat','vs');
fprintf('saved\n');
