function Energy = lc_energy_2d_polar(matrix, K11, K22, K33, deltax, deltay)

Np = matrix;

[dNpdy, dNpdx] = gradient(Np, deltaz, deltay)

Energy = (1/2).*K11.*(dNpdy.*cos(Np)+(-1).*dNpdx.*sin(Np)).^2+(1/2).*K33.*( ...
  dNpdx.*cos(Np)+dNpdy.*sin(Np)).^2;
