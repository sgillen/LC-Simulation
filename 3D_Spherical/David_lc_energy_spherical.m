function [ energy ] = David_lc_energy_spherical(matrix, K11, K22, K33, K24, deltax, deltay, deltaz, wbcMatrix, Ws)

np = matrix(:,:,:,1);
nt = matrix(:,:,:,2);

[dntdx, dntdy, dntdz] = gradient(nt, deltax, deltay, deltaz);
[d2ntdxdx, d2ntdxdy, d2ntdxdz] = gradient(dntdx,deltax,deltay,deltaz);
[~, d2ntdydy, d2ntdydz] = gradient(dntdy,deltax,deltay,deltaz); % 1 repeat
[~, ~, d2ntdzdz] = gradient(dntdz,deltax,deltay,deltaz); % 2 repeats

[dnpdx, dnpdy, dnpdz] = gradient(np, deltax, deltay, deltaz);
[d2npdxdx, d2npdxdy, d2npdxdz] = gradient(dnpdx,deltax,deltay,deltaz);
[~, d2npdydy, d2npdydz] = gradient(dnpdy,deltax,deltay,deltaz); % 1 repeat
[~, ~, d2npdzdz] = gradient(dnpdz,deltax,deltay,deltaz); % 2 repeats

dntdz = unwrap(dntdz,pi,1);
bulkenergy = sum(sum(sum((1/2).*dntdz.^2.*K33+(-1).*Ws.*sin(nt).^2)));

energy = deltax*deltay*deltaz*bulkenergy + deltax*deltay*sum(sum(sum(wbcMatrix.*((-1).*Ws.*sin(nt).^2))));

return
end