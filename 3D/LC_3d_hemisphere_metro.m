k=1;
R=1;
Nx=20;
Ny=Nx;
Nz=20;
xstep=2*R/(Nx-1);
ystep=2*R/(Ny-1);
zstep=2*R/(Nz-1);

[X Y Z]=meshgrid((1:Nx),(1:Ny),(1:Nz));
hemisphere=((X-(Nx+1)/2).^2./((Nx-1)^2/4)+(Y-(Ny+1)/2).^2./((Ny-1)^2/4)+...
            (Z./Nz).^2)<1;

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

inside=hemisphere-edge;

%flat surface. These points always have theta=pi/2;
edge(:,:,1)=0;
flat=hemisphere(:,:,1);

theta=inside.*zeros(Nx,Ny,Nz);
phi=inside.*zeros(Nx,Ny,Nz);

%set BC's
%radial on edge
phi=phi+hemisphere.*atan((Y-(Ny+1)/2)./(X-(Nx+1)/2));
theta=theta+hemisphere.*acos((Z-(Nz+1)/2)./sqrt((X-(Nx+1)/2).^2+(Y-(Ny+1)/2).^2+(Z-(Nz+1)/2).^2));
theta(flat)=pi/2;

xscaled=xstep*(X-(Nx+1)/2);
yscaled=ystep*(Y-(Ny+1)/2);
zscaled=zstep*(Z-(Nz+1)/2);

tfinal=10;
randx=randi(Nx,3*tfinal,1);
randy=randi(Ny,3*tfinal,1);
randz=randi(Nz,3*tfinal,1);

t=1;
while t<=tfinal
    if ~inside(randx(t),randy(t),randz(t))
        randx(t)=[];
        randy(t)=[];
        randz(t)=[];
    else
        t=t+1;
    end
end

randx=randx(1:tfinal);
randy=randy(1:tfinal);
randz=randz(1:tfinal);

for t=1:tfinal
    %Find our location
    xt=randx(t);
    yt=randy(t);
    zt=randz(t);
    %Compute the dot product in spherical coordinates with nearest
    %neighbors    
    cos(theta(xt,yt,zt))*cos(theta(xt,yt,zt))
    
end

%Find nx,ny,nz. Flip so the divergence doesn't go funny if you are doing
%energy calculations.
nx=sin(theta).*cos(phi).*hemisphere.*((xscaled>0)-0.5)*2;
ny=sin(theta).*sin(phi).*hemisphere.*((xscaled>0)-0.5)*2;
nz=cos(theta).*hemisphere;

%The divergence calculation gets goofy on the edges
%For radial configuration, the divergence on the boundaries goes to 2/R,
%so just set those values manually.
% div=hemisphere.*divergence(xscaled,yscaled,zscaled,nx,ny,nz);
% div=div.*(div>0)+2/R.*(div<0);

% divinside=inside.*(((circshift(nx,[0 -1 0])-circshift(nx,[0 1 0]))/(2*xstep)) ...
%     +((circshift(ny,[-1 0 0])-circshift(ny,[1 0 0]))/(2*ystep)) ...
%     +((circshift(nz,[0 0 -1])-circshift(nz,[0 0 1]))/(2*zstep)));

figure(2)
clf
zlook=[1 round(Nz/3) round(2*Nz/3) Nz-3];
ha=tight_subplot(length(zlook),1,0.005,0.005,0.005);
for ii=1:length(zlook)
    axes(ha(ii));
    axis xy;axis square;
    quiver(nx(:,:,zlook(ii)),ny(:,:,zlook(ii)),0)
    axis([0 Nx+1 0 Ny+1])
end

figure(3)
clf
zlook=[1 round(Nz/3) round(2*Nz/3) Nz-3];
ha=tight_subplot(length(zlook),1,0.005,0.005,0.005);
for ii=1:length(zlook)
    axes(ha(ii));
    axis xy;axis square;
    imagesc(sin(2*atan(ny(:,:,zlook(ii))./nx(:,:,zlook(ii)))).^2)
    axis([0 Nx+1 0 Ny+1])
    colormap bone
end

% %Splay energy
% 8*pi*k*R
% splay=k/2*sum(sum(sum(div.^2)))*xstep*ystep*zstep;
% 
% %Bend energy
% [cx cy cz]=curl(xscaled,yscaled,zscaled,nx,ny,nz);
% cx(isnan(cx))=0;
% cy(isnan(cy))=0;
% cz(isnan(cz))=0;
% 
% bend=k/2*sum(sum(sum((cx.*nx+cy.*ny+cz.*nz).^2)))*xstep*ystep*zstep;
% 
% %Twist energy
% ncross=zeros(Nx,Ny,Nz,3);
% ncross(:,:,:,1)=nx;
% ncross(:,:,:,2)=ny;
% ncross(:,:,:,3)=nz;
% 
% ccross=zeros(Nx,Ny,Nz,3);
% ccross(:,:,:,1)=cx;
% ccross(:,:,:,2)=cy;
% ccross(:,:,:,3)=cz;
% 
% crossed=cross(ncross,ccross);
% crossedsquared=inside.*(crossed(:,:,:,1).^2+crossed(:,:,:,2).^2+crossed(:,:,:,3).^2);
% twist=k/2*sum(sum(sum(crossedsquared)))*xstep*ystep*zstep;
% 
% [splay bend twist]