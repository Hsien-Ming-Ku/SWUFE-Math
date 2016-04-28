function [A,b] = nonsymmetric_BBC(n,nu)
% This function is to generate the non-Hermitian linear systems A * x = b
% arising from the BBC Problem proposed by Prof. Bai in [1,2].
% Reference
% 1. Z.-Z. Bai, M. Benzi, F. Chen, Modified HSS iteration methods for a 
%    class of complex symmetric linear systems, Computing, vol. 87, 2010, 
%    pp. 93–111.
% 2. Z.-Z. Bai, On preconditioned iteration methods for complex linear 
%    systems, J. Eng. Math., vol. 93. 2015, pp. 41-60.
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
% Date: 28-4-2016, 14:47 (The Netherlands)
% ---------------------------------------------------------
theta = nu/(2*(n+1));
r1 = [2 -1+theta zeros(1,n-2)];
r2 = [2 -1-theta zeros(1,n-2)];
L = sparse(toeplitz(r2,r1));
I = speye(n,n);
W = kron(L,I) + kron(I,L);
e1 = [1; zeros(n-1,1)]; en = [zeros(n-1,1); 1];
Lc = L - e1*en' - en*e1';
% -- In fact, we think that parameters (i.e, 10 and 9) can be modified by
%    your own choice !!!
T = 10*(kron(Lc,I) + kron(I,Lc)) + 9*kron(e1*en' + en*e1',I); 
A = W + 1i*T;
b = (1 + 1i)*(A*randn(n^2,1));
