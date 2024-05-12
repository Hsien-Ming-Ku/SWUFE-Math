function [A,b] = CDR_threedim(nx,ny,nz,beta,r,opt)
% Cheung W M, Ng M K. Block-circulant preconditioners for systems arising
% from discretization of the three-dimensional convection-diffusion equation[J].
% Journal of Computational and Applied Mathematics, 2002, 140(1):143-158.
hx = 1 / (nx + 1);
hy = 1 / (ny + 1);
hz = 1 / (nz + 1);
Ix = speye(nx,nx);Iy = speye(ny,ny);Iz = speye(nz,nz);
In = speye(nx * ny * nz,nx * ny * nz);
ex = ones(nx,1);ey = ones(ny,1);ez = ones(nz,1);
if opt == 1%中心差分
    a = 6;
    b = -1 - beta(1)*hx / 2;c = -1 + beta(1)*hx / 2;
    d = -1 - beta(2)*hy / 2;e = -1 + beta(2)*hy / 2;
    f = -1 - beta(3)*hz / 2;g = -1 + beta(3)*hz / 2;
    Tx = spdiags([b*ex a*ex c*ex],[-1 0 1],nx,nx)/hx^2;
    Ty = spdiags([d*ey 0*ey e*ey],[-1 0 1],ny,ny)/hy^2;
    Tz = spdiags([f*ez 0*ez g*ez],[-1 0 1],nz,nz)/hz^2;
    A = kron(kron(Tz,Iy),Ix)+kron(kron(Iz,Ty),Ix)+...
        kron(kron(Iz,Iy),Tx)-r*In;
else%迎风格式
    a = 2 + beta(1)*hx + beta(2)*hy + beta(3)*hz;
    b = -1 - beta(1)*hx;c = -1;
    d = -1 - beta(2)*hy;e = -1;
    f = -1 - beta(3)*hz;g = -1;
    Tx = spdiags([b*ex a*ex c*ex],[-1 0 1],nx,nx)/hx^2;
    Ty = spdiags([d*ey 0*ey e*ey],[-1 0 1],ny,ny)/hy^2;
    Tz = spdiags([f*ez 0*ez g*ez],[-1 0 1],nz,nz)/hz^2;
    A = kron(kron(Tz,Iy),Ix)+kron(kron(Iz,Ty),Ix)+...
        kron(kron(Iz,Iy),Tx)-r*In;    
end
 nl = size(A,2);
 x = rand(nl,1);
 b = A*x;
 
