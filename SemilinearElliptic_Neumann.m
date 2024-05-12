function [A,b] = SemilinearElliptic_Neumann(nu,omega,d)
% Hwang T M, Wang W. Analyzing and visualizing a discretized semilinear
% Elliptic problem with Neumann boundary conditions[J]. Numerical Methods
% for Partial Differential Equations, 2002, 18(3):261-279.
delr = 2/(2*nu + 1);
deltheta = 2*pi/omega;
r1 = ones(omega,1);   
Phi = spdiags([-r1 2*r1 -r1],[-1 0 1],omega,omega);
Phi(omega,1) = -1;Phi(1,omega)=-1;
delta = -2 - delr^2/d;
I = speye(omega);
r = zeros(nu,1);
mu = r; beta = r;
for k = 1:nu
    r(k) = (k -1/2)*delr;
    mu(k) = 1/(2*k-1);
    beta(k) = 1/(((k-1/2)^2)*deltheta^2);
end
eta = -1 + mu(nu) - delr^2/d;
a = [delta*ones(nu - 1,1);eta];
b = [0;1 + mu(1:nu-1)]; c = [1 - mu(2:end);0];
Y = spdiags([c a b],[-1 0 1],nu,nu);
A1 = kron(Y,I);
A2 = kron(spdiags(beta,0,nu,nu),Phi);
A = A1 - A2;

nl = size(A,2);
x = rand(nl,1);
b = A*x;
