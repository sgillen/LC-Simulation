function [Elp] = EL_terms_2D_polar(matrix, K11, K22, K33, deltax, deltay, wbcMatrix)


Np = matrix(:,:);
%{
dNp = wrap_gradient(Np);
dNpdx = dNp(1,:);
dNpdy = dNp(2,:);

d2Npdx = wrap_gradient(dNpdx);
d2Npdxdx = d2Npdx(1,:);
d2Npdxdy = dNpdx(2,:);

d2Npdy = wrap_gradient(dNpdy);
d2Npdydx = d2Npdy(1,:);
d2Npdydy = dNpdy(2,:);
%}



Np = matrix(:,:);
[dNpdx, dNpdy] = gradient(Np, deltax, deltay);
[d2Npdxdx, d2Npdxdy] = gradient(dNpdx, deltax, deltay);
[d2Npdydx, d2Npdydy] = gradient(dNpdy, deltax, deltay);


%{
Np = matrix(:,:);

%Take the step size out of gradient and it will simply difference the list.
[dNpdx, dNpdy] = gradient(Np);


%Correct the differences such that the range always falls between
%-pi/2 and pi/2
dNpdx = rem(dNpdx,pi);
dNpdy = rem(dNpdy,pi);

[d2Npdxdx, d2Npdxdy] = gradient(dNpdx);
[~, d2Npdydy] = gradient(dNpdy);


%this was done out of desperation, it did not work
%d2Npdxdx = mod(d2Npdxdx+pi/2,pi) - pi/2;
%d2Npdxdy = mod(d2Npdxdy+pi/2,pi) - pi/2;
%d2Npdydy = mod(d2Npdydy+pi/2,pi) - pi/2;

%Now put the step sizes back in.
dNpdx = dNpdx/deltax;
dNpdy = dNpdy/deltay;

d2Npdxdx = d2Npdxdx/(deltax^2);
d2Npdxdy = d2Npdxdy/deltax/deltay;
d2Npdydy = d2Npdydy/(deltay^2);
%}
Elp = zeros(size(matrix));

Elpb = (1/2).*K33.*((-2).*d2Npdxdx.*cos(Np).^2+(-2).*dNpdx.*dNpdy.*cos(2.*Np)+ ...
  (-2).*d2Npdydy.*sin(Np).^2+(-2).*d2Npdxdy.*sin(2.*Np)+dNpdx.^2.*sin(2.* ...
  Np)+(-1).*dNpdy.^2.*sin(2.*Np))+(1/2).*K11.*(2.*dNpdx.*dNpdy.*cos(2.*Np) ...
  +(-2).*(d2Npdydy.*cos(Np).^2+sin(Np).*((-2).*d2Npdxdy.*cos(Np)+ ...
  dNpdx.^2.*cos(Np)+d2Npdxdx.*sin(Np)))+dNpdy.^2.*sin(2.*Np));


%We can speed things up a bit by not copying into a and Ws and just using
%the bc matrix directly. 
a = wbcMatrix(:,:,1);
Ws = wbcMatrix(:,:,2);
 

Elps = (1/2).*K33.*((-1).*d2Npdxdx+(-1).*(d2Npdxdx+2.*dNpdx.*dNpdy).*cos(2.* ...
  Np)+(-2).*d2Npdydy.*sin(Np).^2+((-2).*d2Npdxdy+dNpdx.^2+(-1).*dNpdy.^2) ...
  .*sin(2.*Np))+(1/2).*K11.*((-1).*d2Npdxdx+(-2).*d2Npdydy.*cos(Np).^2+( ...
  d2Npdxdx+2.*dNpdx.*dNpdy).*cos(2.*Np)+(2.*d2Npdxdy+(-1).*dNpdx.^2+ ...
  dNpdy.^2).*sin(2.*Np))+(-1).*Ws.*sin(2.*((-1).*Np+a));

screen = (~(Ws==0));
Elp(screen) = Elps(screen);
Elp(~screen) = Elpb(~screen);


%Single Constant Approximation 
%Elp = -K11*(d2Npdxdx + d2Npdydy);


end
 