function [Elp] = EL_terms_2D_polar(matrix, K11, K22, K33, deltax, deltay)

Np = matrix(:,:);
[dNpdx, dNpdy] = gradient(Np, deltax, deltay);
[d2Npdxdx, d2Npdxdy] = gradient(dNpdx, deltax, deltay);
[d2Npdydx, d2Npdydy] = gradient(dNpdy, deltax, deltay);

Elp = (1/2).*K33.*((-2).*d2Npdxdx.*cos(Np).^2+(-2).*dNpdx.*dNpdy.*cos(2.*Np)+ ...
  (-2).*d2Npdydy.*sin(Np).^2+(-2).*d2Npdxdy.*sin(2.*Np)+dNpdx.^2.*sin(2.* ...
  Np)+(-1).*dNpdy.^2.*sin(2.*Np))+(1/2).*K11.*(2.*dNpdx.*dNpdy.*cos(2.*Np) ...
  +(-2).*(d2Npdydy.*cos(Np).^2+sin(Np).*((-2).*d2Npdxdy.*cos(Np)+ ...
  dNpdx.^2.*cos(Np)+d2Npdxdx.*sin(Np)))+dNpdy.^2.*sin(2.*Np));

end
 