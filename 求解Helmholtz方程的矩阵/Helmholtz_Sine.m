function [A,P] = Helmholtz_Sine(m,k,p,q)
% Le J, Jin H, Lv X G, Cheng Q S. A preconditioned method for the solution of
% the Robbins problem for the Helmholtz equation[J]. ANZIAM Journal, 2010,
% 52(01):87-100.
h = 1/m;
Im = speye(m);
%求B
r1 = ones(m,1);
B = spdiags([-r1 (4-k^2*h^2)*r1 -r1],[-1 0 1],m,m);
B(m,m-1) = -2;B(m,m) = 4-k^2*h^2 - 2*h*q;
%求C
C = spdiags([-r1 0*r1 -r1],[-1 0 1],m,m);
C(m,m-1) = -2;C(m,m) = - 2*h*p;
%求A
A = kron(Im,B)+kron(C,Im);

%求正弦变换的预条件子P
hatB = spdiags([-r1 (4-k^2*h^2)*r1 -r1],[-1 0 1],m,m);
hatB(m,m-1) = -2;
P = kron(Im,hatB)+kron(C,Im);
