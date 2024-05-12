function [A,b]=semi_circulant(m1,m2,omega,Sigma)
% Hemmingsson L. A semi-circulant preconditioner for the convection-diffusion
% equation[J]. Numerische Mathematik, 1998, 81(2):211-248.
h1=1/(m1+1);h2=1/(m2+2);
alpha1=omega(1)/h1;alpha2=omega(2)/h2;
beta1=2*Sigma/(h1^2);beta2=2*Sigma/(h2^2);
%求A1
r1 = ones(m1,1);
A1 = spdiags([(-alpha1-beta1)*r1 2*beta1*r1 (alpha1-beta1)*r1],...
    [-1 0 1],m1,m1);
%求A2
r2 = ones(m2,1);
A2 = spdiags([(-alpha2-beta2)*r2 2*beta2*r2 (alpha2-beta2)*r2],...
    [-1 0 1],m2,m2);
%求A
Im1 = speye(m1,m1);Im2 = speye(m2,m2);
A = kron(Im2,A1)+kron(A2,Im1);
%求M,M是一个预处理器,M^(-1)*A*x = M^(-1)*b
C1 = A1;C1(m1,1)=alpha1-beta1;C1(1,m1)=-alpha1-beta1;
M = kron(Im2,C1)+kron(A2,Im1);

nl = size(M,2);
x = rand(nl,1);

b = A*x;
