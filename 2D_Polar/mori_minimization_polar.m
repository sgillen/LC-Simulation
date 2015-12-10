
%function [phi] = mori_minimization_polar(Ws)
%% initalize variables we'll need
%clear all
% how many steps to take?
numsteps = 10000;
% save this number of orientations to capture dynamics
frames = 100;
% how large of mesh
gridsize = [10 10];
% how long should each step represent for dynamics? Measured in seconds?
timestep = 0.000001;
% the LC molecules realign at a rate determined by rotational viscosity.
gamma = 0.08; % rotational viscosity of 80 centiPoise (MBBA at room temp)
% how large is the sample area?
cellsize = [2e-6, 2e-6];
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


%% begin with some initial configuration

% %  randomly oriented
nMatrix = randn(gridsize(1),gridsize(2));

% % fixed orientation
% nMatrix = ones(gridsize(1),gridsize(2),2);
nMatrix(:,:) = 0;

%% set some boundary conditions. these should be enforced at each step.

%create a matrix representing the weak BC's enforced by the modifed L terms
wbcMatrix = create_weak_bcs(size(nMatrix));

%enforce the strong boundaries
[nMatrix] = enforce_2d_polar_sbcs(nMatrix);

% Sort of like normalizing, but smarter about it, I'm not entirley sure
% this is actually working.. 

nMatrix = unwrap(nMatrix,pi,1);

% initialize an array to store snapshots of the nMatrix in time.
nMnMatrixSnapshot = repmat(zeros(size(nMatrix)),1,1,frames);

% set the first snapshot to the inital condition
nMatrixSnapshot(:,:,1) = nMatrix;

% take a picture whenever we hit one of the following time periods
snapshot = round(linspace(1,numsteps,frames))';
snapnum = 2;

avgEnergy = zeros(1,frames);
avgEnergy(1) = mean2(lc_energy_2d_polar(nMatrix,K11,K22,K33,dx,dy,wbcMatrix));


for ii = 2:numsteps
    %calulate step sizes
    [ELp] = el_terms_2d_polar(nMatrix,K11,K22,K33,dx,dy,wbcMatrix);
    nMatrix = nMatrix - timestep./gamma.*(ELp);
    
    %enforce BCs
    nMatrix = enforce_2d_polar_sbcs(nMatrix);
    
    % renormalize 
    nMatrix = unwrap(nMatrix,pi);
    %nMatrix = mod(nMatrix,2*pi);
    
    
    if any(ii==snapshot)
        fprintf('%d%%\n',round(100*ii/numsteps))
        fprintf('%d\n',timestep/gamma*max(max(max(ELp))))
        %avg(snapnum) = mean2(lc_energy_2D_Cartesian(nMatrix, K11, K22, K33, dx,dy)) ;
        avgEnergy(snapnum) = mean2(lc_energy_2d_polar(nMatrix,K11,K22,K33,dx,dy,wbcMatrix));
        fprintf('energy: %d\n', avgEnergy(snapnum)) 
        nMatrixSnapshot(:,:,snapnum) = nMatrix;
        snapnum = snapnum + 1;
    end
end

figure(1)
clf
q2 = quiver(cos(nMatrix), sin(nMatrix)) ;
set(q2,'ShowArrowHead','on')

figure(2)
plot(avgEnergy(2:end))
xlabel('frame')
ylabel('average energy')

%% Label and save everything

%let's contruct some strings so we can save things
home_dir = '/home/gillen/';
save_dir = '/Documents/Computation/saved_outputs/2D_unwrap_test/';
%name_dir =  sprintf('Time_step:%d__frames__WS:%d', wbcMatrix(1,end),timestep,frames);
%name_dir = sprintf('Time_Step:%d__Date:%s/' ,timestep, datestr(datetime('now')));
name_dir  = sprintf('0start/');


eng_fig_name = 'energy.jpg';
director_fig_name = 'director.jpg';
final_phi = 'phi';

%not sure why but matlab tries to save to my home directory instead of the
%location of the script. specifying the absolute path fixes this, but you'll
%need to change it if you run this on some other system. 

%if(~exists('strcat(home_dir,save_dir,name_dir)'))
    mkdir(strcat(home_dir,save_dir,name_dir));
%end
save(strcat(home_dir,save_dir,name_dir,'nMatrix.mat'),'nMatrixSnapshot');
save('n2d.mat' , 'nMatrixSnapshot');
save(strcat(home_dir,save_dir,name_dir,'avgEnergy'), 'avgEnergy');
save(strcat(home_dir,save_dir,name_dir,'workspace'))%without arguments save saves the whole workspace, This shouldn't bee too h
saveas(figure(2), strcat(home_dir,save_dir,name_dir,eng_fig_name));
saveas(figure(1), strcat(home_dir,save_dir,name_dir,director_fig_name));

phi = nMatrix(1,end);
save(strcat(home_dir,save_dir,name_dir,'final_phi'), 'phi');

%end
