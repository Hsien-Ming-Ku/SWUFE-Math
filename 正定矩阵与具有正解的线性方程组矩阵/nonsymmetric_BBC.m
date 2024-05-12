function [A,b1] = nonsymmetric_BBC(n,nu)
% This function is to generate the non-Hermitian linear systems A * x = b
% arising from the BBC Problem proposed by Prof. Bai in [1,2].
% Bai Z Z, Benzi M, Chen F. Modified HSS iteration methods for a class of
% complex symmetric linear systems[J]. Springer Vienna, 2010, 87(3):93-111.
% ---------------------------------------------------------
theta = nu/(2*(n+1));
r = ones(n,1);
L = spdiags([(-1-theta)*r 2*r (-1+theta)*r],[-1 0 1],n,n);
I = speye(n,n);
W = kron(L,I) + kron(I,L);
e1 = [1; zeros(n-1,1)]; en = [zeros(n-1,1); 1];
Lc = L - e1*en' - en*e1';
% -- In fact, we think that parameters (i.e, 10 and 9) can be modified by
%    your own choice !!!
T = 10*(kron(Lc,I) + kron(I,Lc)) + 9*kron(e1*en' + en*e1',I); 
A = W + 1i*T;
b1 = (1 + 1i)*(A*ones(n^2,1));
%nl = size(A,2);
%b = ones(nl,1);
