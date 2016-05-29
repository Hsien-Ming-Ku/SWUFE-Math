function [A,b] = CDR_threedim(nx,ny,nz,beta,r,opt)
%
% beta = (tau, simga, mu)  % convection coefficients, when tau. sigma, mu
%                            become larger, the discretized matrix A will
%                            be more ill-conditioned.
%
%
% Reference
% 1. W.M. Cheung, M.K. Ng, Block-circulant preconditioners for systems 
%    arising from discretization of the three-dimensional 
%    convection-diffusion equation, J. Comput. Appl. Math., 140 (2002), 
%    pp. 143-158.
% 2. E. Jiang, Algorithm for solving shifted skew-symmetric linear system,
%    Front. Math. China., 2 (2007), pp. 227-242.
% -----------------------------------------------------------------------
hx = 1/(1 + nx); hy = 1/(1 + ny); hz = 1/(1 + nz);
tau = beta(1); sigma = beta(2); mu = beta(3);
Ix = speye(nx, nx); Iy = speye(ny, ny); Iz = speye(nz, nz);
In = speye(nx*ny*nz,nx*ny*nz);
ex = ones(nx,1); ey = ones(ny,1); ez = ones(nz,1); 
if opt == 1 % centered differences for the first-order derivatives
   a = 2;  b = -1 - tau*hx/2; d = -1 - sigma*hy/2; f = -1 - mu*hz/2;
   c = -1 + tau*hx/2; e = -1 + sigma*hy/2;  g = -1 + mu*hz/2;
   Tx =  spdiags([b*ex a*ex c*ex],[-1 0 1],nx,nx)/hx^2;
   Ty =  spdiags([d*ey a*ey e*ey],[-1 0 1],ny,ny)/hy^2;
   Tz =  spdiags([f*ez a*ez g*ez],[-1 0 1],nz,nz)/hz^2;
   A = kron(kron(Tz,Iy),Ix) + kron(kron(Iz,Ty),Ix) + ...
       kron(kron(Iz,Iy),Tx) - r*In;
   x = linspace(hx,1 - hx,nx); y = linspace(hy,1 - hy,ny); 
   z = linspace(hz,1 - hz,nz);
   sol = kron(kron(x.*(1 - x),y.*(1 - y)),z.*(1 - z))';
   b = A*sol;
else % backward difference approximations for the first-order derivatives
   ax = 2 + tau*hx;  ay = 2 + sigma*hy; az = 2 + mu*hz;
   b = -1 - tau*hx; d = -1 - sigma*hy; 
   f = -1 - mu*hz;  c = -1;  e = -1;   g = -1;
   Tx = spdiags([b*ex ax*ex c*ex],[-1 0 1],nx,nx)/hx^2;
   Ty = spdiags([d*ey ay*ey e*ey],[-1 0 1],ny,ny)/hy^2;
   Tz = spdiags([f*ez az*ez g*ez],[-1 0 1],nz,nz)/hz^2;
   A = kron(kron(Tz,Iy),Ix) + kron(kron(Iz,Ty),Ix) + ...
       kron(kron(Iz,Iy),Tx) - r*In;
   x = linspace(hx,1 - hx,nx); y = linspace(hy,1 - hy,ny); 
   z = linspace(hz,1 - hz,nz);
   sol = kron(kron(x.*(1 - x),y.*(1 - y)),z.*(1 - z))';
   b = A*sol;
end
