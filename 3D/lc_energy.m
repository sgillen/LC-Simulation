%% This function should return the energy for the 3D Cartesian Oseen Frank Rep 

function energy = lc_energy(matrix, K11, K22, K33, deltax, deltay, deltaz)

Nx = matrix(:,:,:,1);
Ny = matrix(:,:,:,2);
Nz = matrix(:,:,:,3);

[dNxdx, dNxdy, dNxdz] = gradient(Nx, deltax, deltay, deltaz);
[dNydx, dNydy, dNydz] = gradient(Ny, deltax, deltay, deltaz);
[dNzdx, dNzdy, dNzdz] = gradient(Nz, deltax, deltay, deltaz);


energy = (1/2).*(dNxdx+dNydy+dNzdz).^2.*K11+(1/2).*K22.*(((-1).*dNydz+ ...
dNzdy).*Nx+(dNxdz+(-1).*dNzdx).*Ny+((-1).*dNxdy+dNydx).*Nz).^2+( ...
1/2).*K33.*(((dNxdz+(-1).*dNzdx).*Nx+(dNydz+(-1).*dNzdy).*Ny).^2+( ...
(dNxdy+(-1).*dNydx).*Ny+(dNxdz+(-1).*dNzdx).*Nz).^2+(((-1).*dNxdy+ ...
dNydx).*Nx+(dNydz+(-1).*dNzdy).*Nz).^2);

end