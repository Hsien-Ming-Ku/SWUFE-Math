function [A, b] = SemilinearElliptic_Neumann(nu,omega,d)
% Reference
% 1. T.-M. Hwang and W. Wang, Analyzing and visualizing a discretized 
%    semilinear elliptic problem with Neumann boundary conditions, Numer. 
%    Methods for Partial Differential Equations, Vol. 18, pp 261-279, 2002.
% ---------------------------------------------------------
% Copyright by 
% Mr. Xian-Ming Gu, who works as a PhD student in School of Mathematical
% Sciences, University of Electronic Science and Technology of China 
% (UESTC). From 1 August, 2014 to 1 Augest, 2016, He also work in the 
% Johann Bernoulli Institute for Mathematics and Computer Science (JBI), 
% University of Groningen (RUG). Prof. Ting-Zhu Huang (UESTC) and Dr. Bruno
% Carpentieri (RUG) are Xian-Ming's supervisor for PhD program.
% --------------------------------------------------------------------
% Email: guxianming@live.cn
% URL: https://github.com/Hsien-Ming-Ku
% Date: 28-4-2016, 13:30 (The Netherlands)
% ---------------------------------------------------------
delr = 2/(2*nu + 1);
deltheta = 2*pi/omega;
r1 = [2 -1 zeros(1,omega - 3) -1];    r2 = [2 -1 zeros(1,omega - 3) -1];
Phi = sparse(toeplitz(r2,r1));
delta = -2 - delr^2/d;
I = speye(omega);
r = zeros(1,nu);
mu = r; beta = r;
for k = 1:nu
    r(k) = (k -1/2)*delr;
    mu(k) = 1/(2*k-1);
    beta(k) = 1/(((k-1/2)^2)*deltheta^2);
end
eta = -1 + mu(nu) - deltheta^2;
Y = diag([delta*ones(1,nu - 1) eta]) + diag(1 + mu(1:nu-1),1) + ...
    + diag(1 - mu(2:end),-1);
A1 = kron(Y,I);
A2 = kron(-diag(beta),Phi);
A = A1 + A2;
nl = size(A,2);
b = A*ones(nl,1);
