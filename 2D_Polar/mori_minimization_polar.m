%% initalize variables we'll need
clear all
% how many steps to take?
numsteps = 40000;
% save this number of orientations to capture dynamics
frames = 50;
% how large of mesh
gridsize = [10 10];
% how long should each step represent for dynamics? Measured in seconds.
timestep = 0.00025;
% the LC molecules realign at a rate determined by rotational viscosity.
gamma = 0.08; % rotational viscosity of 80 centiPoise (MBBA at room temp)
% how large is the sample area?
cellsize = [20e-6, 20e-6];
% calculate size of each mesh point
dx = cellsize(1)/(gridsize(1)-1); % minus one because matlab starts indexing at 1
dy = cellsize(2)/(gridsize(2)-1);
% elastic coefficients for LC, in SI units of N*m.
K11 = 12e-12; % 5CB
K22 = 5e-12;
K33 = 12e-12;
%% define an XY grid for plotting
%Always use ndgrid and never meshgrid meshgrid will flip the X and Y coordinates.
[X,Y] = ndgrid(1:gridsize(1),1:gridsize(2));

% Create a matrix to represent the director in our grid. For each grid
% point, there is an nx and an ny value. Although the n-vector is
% normalized (which means there is only one free variable) and we only need
% the angular orientation of the n-vector, you can get numerical
% instabilities from the trig functions. It's easier (but a little more
% memory-intensive) to just calculate nx and ny separately and normalize at
% each step. For 2D, it may well be best to just calculate the angle. But
% the method we use here generalizes to 3D more easily as well.



%% begin with some initial configuration

% %  randomly oriented
nMatrix = randn(gridsize(1),gridsize(2));

% % fixed orientation
% nMatrix = ones(gridsize(1),gridsize(2),2);
% nMatrix(:,:,2) = 1;

%% set some boundary conditions. these should be enforced at each step.

nMatrix = enforce_2d_polar_bcs(nMatrix);

% normalize. this should also be done at each step.
% I'm not sure why but the way I'm doing things doesn't like angles outside
% the range 0 < phi < 2pi ,this could be a sign that something else is
% wrong though. 
nMatrix = mod(nMatrix,2*pi);


% initialize an array to store snapshots of the nMatrix in time.
nMnMatrixSnapshot = repmat(zeros(size(nMatrix)),1,1,frames);

% set the first snapshot to the inital condition
nMatrixSnapshot(:,:,1) = nMatrix;

% take a picture whenever we hit one of the following time periods
snapshot = round(linspace(1,numsteps,frames))';
snapnum = 2;

[ELp] = el_terms_2d_polar(nMatrix,K11,K22,K33,dx,dy);

avgEnergy = zeros(1,frames);
avgEnergy(1) = mean2(lc_energy_2d_polar(nMatrix,K11,K22,K33,dx,dy));

for ii = 2:numsteps
    % calulate step sizes
    [ELp] = el_terms_2d_polar(nMatrix,K11,K22,K33,dx,dy);
    nMatrix = nMatrix - timestep./gamma.*ELp;
    
    %enforce BCs
    nMatrix = enforce_2d_polar_bcs(nMatrix);
    
    % renormalize 
    nMatrix = mod(nMatrix,2*pi);
    %nMatrix = nMatrix./repmat(sqrt(sum(nMatrix.^2,3)),[1 1 2]);
    
    if any(ii==snapshot)
        fprintf('%d%%\n',round(100*ii/numsteps))
        fprintf('%d\n',timestep/gamma*max(max(max(ELp))))
        %avg(snapnum) = mean2(lc_energy_2D_Cartesian(nMatrix, K11, K22, K33, dx,dy)) ;
        avgEnergy(snapnum) = mean2(lc_energy_2d_polar(nMatrix,K11,K22,K33,dx,dy));
        fprintf('energy: %d\n', avgEnergy(snapnum)) 
        nMatrixSnapshot(:,:,snapnum) = nMatrix;
        snapnum = snapnum + 1;
    end
end

figure(1)
clf
q2 = quiver(sin(nMatrix), cos(nMatrix)) ;
set(q2,'ShowArrowHead','on')

figure(2)
plot(avgEnergy(2:end))
xlabel('frame')
ylabel('average energy')

%% Label and save everything

%let's contruct some strings so we can save things
home_dir = '/home/gillen/';
save_dir = '/Documents/Computation/saved_outputs/OF_2D_Linear_twist_polar/';
name_dir = sprintf('Twist:%d-%d/Time_step:%d__frames:%d__Date:%s/', nMatrix(1,1)*180/pi, nMatrix(1,end)*180/pi, timestep, frames, datestr(datetime('now')));


%not sure why but matlab tries to save to my home directory instead of the
%location of the script. specifying the absolute path fixes this, but you'll
%need to change it if you run this on some other system. 

%if(~exists('strcat(home_dir,save_dir,name_dir)'))
    mkdir(strcat(home_dir,save_dir,name_dir));
%end
save(strcat(home_dir,save_dir,name_dir,'nMatrix.mat'),'nMatrixSnapshot');
save(strcat(home_dir,save_dir,name_dir,'aveEnergy'), 'avgEnergy');
save(strcat(home_dir,save_dir,name_dir,'workspace')) %without arguments save saves the whole workspace, This shouldn't bee too h