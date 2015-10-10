%% initalize variables we'll need
% how many steps to take?
numsteps = 40000;
% save this number of orientations to capture dynamics
frames = 50;
% how large of mesh?
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
nMatrix = randn(gridsize(1),gridsize(2),2);

% % fixed orientation
% nMatrix = ones(gridsize(1),gridsize(2),2);
% nMatrix(:,:,2) = 1;

%% set some boundary conditions. these should be enforced at each step.
% enforce BC's
nMatrix(:,1,2) = 0;
nMatrix(:,1,1) = 1;
    
nMatrix(:,end,2) = .5;
nMatrix(:,end,1) = .5;




% nMatrix(:,1,2) = 0;
% nMatrix(:,end,2) = 0;
% nMatrix(1,1,1) = 1; nMatrix(1,1,2) = 1;
% nMatrix(end,1,1) = 1; nMatrix(end,1,2) = -1;
% nMatrix(1,end,1) = 1; nMatrix(1,end,2) = -1;
% nMatrix(end,end,1) = 1; nMatrix(end,end,2) = 1;

% normalize. this should also be done at each step.
nMatrix = nMatrix./repmat(sqrt(sum(nMatrix.^2,3)),[1 1 2]);

% initialize an array to store snapshots of the nMatrix in time.
nMatrixSnapshot = repmat(zeros(size(nMatrix)),1,1,1,frames);

% set the first snapshot to the inital condition
nMatrixSnapshot(:,:,:,1) = nMatrix;

% take a picture whenever we hit one of the following time periods
snapshot = round(linspace(1,numsteps,frames))';
snapnum = 2;

[ELx, ELy] = EL_terms_2D(nMatrix,K11,K22,K33,dx,dy);   

for ii = 2:numsteps
    % calulate step sizes
    [ELx, ELy] = EL_terms_2D(nMatrix,K11,K22,K33,dx,dy);
    nMatrix(:,:,1) = nMatrix(:,:,1) - timestep./gamma.*ELx;
    nMatrix(:,:,2) = nMatrix(:,:,2) - timestep./gamma.*ELy;
    
    % enforce BC's
    nMatrix(:,1,2) = 0;
    nMatrix(:,1,1) = 1;
    
    nMatrix(:,end,2) = .5;
    nMatrix(:,end,1) = .5;
    
    % renormalize
    nMatrix = nMatrix./repmat(sqrt(sum(nMatrix.^2,3)),[1 1 2]);
    
    if any(ii==snapshot)
        fprintf('%d%%\n',round(100*ii/numsteps))
        fprintf('%d\n',timestep/gamma*max(max(max(ELx))))
        avg(snapnum) = mean2(lc_energy_2D_Cartesian(nMatrix, K11, K22, K33, dx,dy)) ;
        nMatrixSnapshot(:,:,:,snapnum) = nMatrix;
        snapnum = snapnum + 1;
    end
end

figure(1)
clf
q2 = quiver(X-nMatrix(:,:,1)./2,Y-nMatrix(:,:,2)./2,nMatrix(:,:,1),nMatrix(:,:,2));
set(q2,'ShowArrowHead','on')

figure(2)
plot(avg(3:end))
xlabel('frame')
ylabel('average energy')

%% Label and save everything

%let's contruct some strings so we can save things
home_dir = '/home/gillen/';
save_dir = '/Documents/Computation/saved_outputs/OF_2D_Linear_twist/';
name_dir = sprintf('Twist:%d-%d/Date:%s/', atand(nMatrix(1,1,2)/nMatrix(1,1,1)), atand(nMatrix(1,end,2)/nMatrix(1,end,2)),datestr(datetime('now')));


%not sure why but matlab tries to save to my home directory instead of the
%location of the script. specifying the absolute path fixes this, but you'll
%need to change it if you run this on some other system. 

%if(~exists('strcat(home_dir,save_dir,name_dir)'))
    mkdir(strcat(home_dir,save_dir,name_dir));
%end
save(strcat(home_dir,save_dir,name_dir,'nMatrix.mat'),'nMatrixSnapshot');
save(strcat(home_dir,save_dir,name_dir,'aveEnergy'), 'avg');
save(strcat(home_dir,save_dir,name_dir,'workspace')) %without arguments save saves the whole workspace