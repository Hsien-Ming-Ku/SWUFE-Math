function [A,P] = Helmholtz_Fourier(N,a,k,alpha,beta)
% Gryazin Y. Preconditioned Krylov subspace methods for sixth order compact
% approximations of the Helmholtz equation[J]. ISRN Computational 
% Mathematics, 2012, 2014:1-15.
h = a/(N+1);
I = speye(N,N);

l = zeros(N,1);
e = ones(N,1);

%求预条件子Ap
lambeta1 = [ones(N-1,1);1+beta];
lambeta2 = [1+beta;ones(N-1,1)];
LamBN = spdiags([lambeta1 l lambeta2],[-1 0 1],N,N);
d1 = [alpha;zeros(N-2,1);alpha];
D = spdiags(d1,0,N,N);
A1 = LamBN-2*I;A2 = LamBN-(2-k^2*h^2)*I+D;
P = kron(kron(A1,I),I)+kron(kron(I,A1),I)+kron(kron(I,I),A2);  
%求A
Lam0N = spdiags([e l e],[-1 0 1],N,N);
Lam00N = Lam0N-2*I;Lam10N = Lam0N-4*I;
A = kron(kron(Lam00N,I),I)+kron(kron(I,Lam00N),I)+kron(kron(I,I),Lam00N)+...
1/6*(1+h^2*k^2/30)*(kron(kron(Lam0N,Lam10N),I)+kron(kron(I,Lam10N),Lam10N)-...
4*kron(kron(I,I),I))+1/30*kron(kron(Lam10N,Lam10N),Lam10N)+...
(h^2*k^2-h^4*k^4/12+h^6*k^6/360)*kron(kron(I,I),I);
