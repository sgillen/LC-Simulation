function Energy = lc_energy_2D_polar(matrix, K11, K22, K33, deltax, deltay)

Np = matrix;

[dNpdx, dNpdy] = gradient(Np, deltax, deltay);

Energy = (1/2).*K11.*(dNpdy.*cos(Np)+(-1).*dNpdx.*sin(Np)).^2+(1/2).*K33.*( ...
  dNpdx.*cos(Np)+dNpdy.*sin(Np)).^2;

end
