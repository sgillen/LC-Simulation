 %% Mori Minimization Polar,

%{

 this program is 3D simulation of nematics liquid crystals, we follow the 
 model set out by mori et al. 

 TODO: link to mori paper. 

 This file is the entry point to the program, press run and off you go.

 To modify strong boundaries look at enforce_3d_spherical_sbcs.m

 To modify weak boundaries look at create_spherical_wbcs.m

 The energy plot is calculated with the expression inside
 lc_energy_spherical.m, but this is not the expression used for
 minimization. 

 anything else you may want to modify (step size, grid size, initial
 condition etc. is going to be in the first two blocks of this script.

 I recommned plotting with plot_3d, you can either plot the whole grid
 "plot_3d(nMatrix)" for example. or only certain slices
 plot_3d(nMatrix(:,:,z,:))
 also worth pointing out is the variable nMatrixFrames which holds all the snapshots
 taken of the matrix. 

%}


%% initialize a grid and other parameters
%how long (in seconds) should each timestep be?
timeStep = 0.0005; 
%number of time steps to take
numSteps = 1000; 
%number of time frames to save
frames = 100;  

%how big is the grid (there is a director at every point)
grid = [2 2 10];
%how big is the total grid (meters?)
cellSize = [2e-6, 2e-6, 2e-6];
%spacing between directors
dx = cellSize(1)/grid(1);
dy = cellSize(2)/grid(2);
dz = cellSize(3)/grid(3);
%Creates a numbered grid thing
[x, y, z] = ndgrid(1:grid(1),1:grid(2),1:grid(3));

%bend, twist, and splay constants
% k11 = 12e-12;
% k22 = 5e-12;
% k33 = 12e-12;
% k24 = 3.1*k11;

k11 = 12e-12;
k22 = 12e-12;
k33 = 12e-12;
k24 = 12e-12;

Ws = 1e-6;
% finalGuess = lastthetas(end) * pi/180;
finalGuess = 45 * pi/180;
finalGuessArray = linspace(pi/4, finalGuess,grid(3));
finalGuessArray = repmat(finalGuessArray(:),[1,grid(1),grid(2)]);
finalGuessArray = shiftdim(finalGuessArray,1);% + 0.002-0.004*randn(grid(1),grid(2),grid(3));
%rotational viscosity of 0.8 Poise
gamma = 0.08; 


%each point in the matrix has an x y and z component, nmatrix(:,:,:,1) is 
%phi, which is the azimuthal angle and nmatrix(:,:,:,2) is theta, the
%elevation angle. 

% nMatrix = zeroes(grid(1),grid(2),grid(3),2);
% Define a logical matrix that is 1 on weak boundaries.
wbcMatrix = zeros(grid(1),grid(2),grid(3));
wbcMatrix(:,:,end) = 1;
% Define a logical matrix that is 0 on strong boundaries.
strMatrix = ones(grid(1),grid(2),grid(3));
strMatrix(:,:,1) = 0;

%% begin with some initial configuration

% %Give an initial guess. 
% nMatrix(:,:,:,1) = 0;
% nMatrix(:,:,:,2) = finalGuessArray;
% %Initialize the strong boundary.
% nMatrix(:,:,1,2) = pi/4;
% initialize arrays that we can use to watch the system evolve
nMatrixFrames = repmat(zeros(size(nMatrix)),1,1,1,1,frames);
nMatrixFrames(:,:,:,:,1) = nMatrix;

% take a picture whenever we hit one of the following time periods
maxStep = zeros(1,frames);
snapShot = round(linspace(1,numSteps,frames))';
snapNum = 1;

avgEnergy = zeros(1,frames);
avgEnergy(1,snapNum) = David_lc_energy_spherical(nMatrix, k11, k22, k33, k24, dx, dy, dz, wbcMatrix,Ws);

snapNum = 2;

%% Do the calculations 
for ii = 2:numSteps
    nMatrix(:,:,:,1) = mod(nMatrix(:,:,:,1),2*pi); 
    nMatrix(:,:,:,2) = mod(nMatrix(:,:,:,2),pi);

    %add noise
    nMatrix(:,:,:,2) = nMatrix(:,:,:,2) + strMatrix.*(0.00001*acos(2*rand(size(nMatrix(:,:,:,2)))-1));

    %grab the EL terms
    [elp, elt]= David_EL_terms(nMatrix,k11,k22,k33,k24,dx,dy,dz,wbcMatrix,Ws);
    nMatrix(:,:,:,1) = nMatrix(:,:,:,1) - strMatrix.*(timeStep/gamma*elp);
    nMatrix(:,:,:,2) = nMatrix(:,:,:,2) - strMatrix.*(timeStep/gamma*elt);    
    
    %save at the certain points we care about. 
    if any(ii==snapShot)      
       avgEnergy(snapNum) = David_lc_energy_spherical(nMatrix, k11, k22, k33,k24, dx, dy, dz,wbcMatrix,Ws);
        fprintf('step %d/%d      avg = %d  final theta = %f\n', snapNum, frames, avgEnergy(snapNum), nMatrix(1,1,end,2)*180/pi);
        nMatrixFrames(:,:,:,:,snapNum) = nMatrix;
        snapNum = snapNum + 1;        
    end
end

%% analysis

%I wrote this one myself, make sure it's in matlabs path.
% plot_3d(nMatrix);

%[vx,vy,vz] = sph2cart(vt,vp,1);
%questionable whether this is really working 
%q3 = quiver3(X-vx(3:end)./2,Y-vy(3:end)./2,Z-vz(3:end)./2,vx,vy,vz);
%set(q3,'ShowArrowHead','off')

figure(2)
clf
plot(avgEnergy(2:end));

lastthetas=nMatrix(1,1,:,2)*180/pi;
lastthetas=lastthetas(:);
figure(3)
clf
plot(lastthetas)
%octant(:,:,:,1) = sqrt(((x-1)/(grid(1)-1)).^2 + ((y-1)/(grid(2)-1)).^2 + ((z-1)/(grid(3)-1)).^2) < 1;
%octant(:,:,:,2) = sqrt(((x-1)/(grid(1)-1)).^2 + ((y-1)/(grid(2)-1)).^2 + ((z-1)/(grid(3)-1)).^2) < 1;
%uMatrix(~octant) = pi/2;

%% let's throw this to mathematica for more analysis
% save('~/n3d.mat','nMatrixFrames');
% fprintf('saved\n');