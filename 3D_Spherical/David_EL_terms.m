function [ ELp, ELt ] = David_EL_terms(matrix, K11, K22, K33, K24, deltax, deltay, deltaz, wbcMatrix, Ws)

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

d2ntdzdz = unwrap(d2ntdzdz,pi,1);

ELt = (-1).*d2ntdzdz.*K33+(-1).*Ws.*sin(2.*nt);

ELt = ELt - Ws*sin(2*nt).*wbcMatrix;% Weak boundary
ELp = zeros(size(ELt));

return
end