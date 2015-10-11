%This function should return the enerfy for the 2D cartesian case

function energy = lc_energy(matrix, K11, K22, K33, deltax, deltay)


Nx = matrix(:,:,1);
Ny = matrix(:,:,2);

[dNxdx, dNxdy] = gradient(Nx, deltax, deltay);
[dNydx, dNydy] = gradient(Ny, deltax, deltay);


energy = (1/2).*(dNxdx+dNydy).^2.*K11+(1/2).*(dNxdy+(-1).*dNydx).^2.*K33.*( ...
Nx.^2+Ny.^2);


end