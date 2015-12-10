function Energy = lc_energy_2D_polar(matrix, K11, K22, K33, deltax, deltay, wbcMatrix)


Np = matrix;
%these copies are totally for convience.
a = wbcMatrix(:,:,1);
Ws = wbcMatrix(:,:,2);


screen = ~(Ws == 0);


[dNpdx, dNpdy] = gradient(Np, deltax, deltay);

%One Constant approx
%Energy = (1/2)*K11*(dNpdx.^2 + dNpdy.^2);

%Two Constant, should be the same as above if K11 = K33

Energyb = (1/2).*K11.*(dNpdy.*cos(Np)+(-1).*dNpdx.*sin(Np)).^2+(1/2).*K33.*( ...
  dNpdx.*cos(Np)+dNpdy.*sin(Np)).^2;

%Energy of the weak boundaries. 


Energys = (1/2).*K11.*(dNpdy.*cos(Np)+(-1).*dNpdx.*sin(Np)).^2+(1/2).*K33.*( ...
  dNpdx.*cos(Np)+dNpdy.*sin(Np)).^2+(-1).*Ws.*(cos(Np).*cos(a)+sin( ...
  Np).*sin(a)).^2;

Energy(screen) = Energys(screen);
Energy(~screen) = Energyb(~screen);


end
