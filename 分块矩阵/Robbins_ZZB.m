function [A,M,L] = Robbins_ZZB(n,q,alpha)
% Bai Z Z, Ng M K. On inexact preconditioners for nonsymmetric matrices[J].
% SIAM Journal on Scientific Computing, 2005, 26(5):1710-1724.
h = 1/n;
r = q*h/2;
r1 = (-1 - r)*ones(n,1); r2 = (-1 + r)*ones(n,1);
Tn = spdiags([r1 2*ones(n,1) r2],[-1 0 1],n,n);
In = speye(n);
en = In(:,n);
x = h:h:1; y = h:h:1;
a = 2^(alpha -1)*sin(y'); c = 2^(alpha -1)*exp(x');
D1 = spdiags(a,0,n,n);  D2 = spdiags(c,0,n,n);
%求A
A = kron(Tn,In) + kron(In,Tn) - 2*h*kron(D1,en*en') - 2*h*kron(en*en',D2);

%求L和M
r_1 = (-1 - r)*ones(n-1,1); r_2 = (-1 + r)*ones(n-1,1);
T_n = spdiags([r_1 2*ones(n-1,1) r_2],[-1 0 1],n-1,n-1);
I_n = speye(n-1);e_n = I_n(:,n-1);
y_1 = h:h:1-h;a_1 = 2^(alpha -1)*sin(y_1'); 
D_1 = spdiags(a_1,0,n-1,n-1);
B = kron(T_n,In)+kron(I_n,Tn)-2*h*kron(D_1,en*en');
F = (-1+r)*kron(e_n,In);E = (-1-r)*kron(e_n',In);
C = 2*In+Tn-2*h*D2-2*h*a(n)*(en*en');
%求L
O = zeros(n*(n-1),n);L = [B O;E C];
%求M
M = [B F;E C];



