function [A,b,x1] = Robbins_ZZB(n,q,alpha)
% This function is used to generate the nonsymmetric linear systems A*x = b
% coming from the finite difference discretization of 2D Robbins Problem in 
% the domain [0,1]\times [0,1], refer to [1,2]. Some notes should be stated 
% here:
% 1) alpha >= 1;
% 2) when q becomes larger, the linear systems becomes more nonsymmetric;
% 3) when alpha becomes large, the linear systems become more indefinite.
% ---------------------------------------------------------------------
% Reference
% 1. Z.-Z. Bai and M.K. Ng, On inexact preconditioners for nonsymmetric
%    matrices, SIAM J. Sci. Comput., vol. 26, no. 5, 2005, pp. 1710-1724.
% 2. W. Pickering, P. Harley, FFT solution of the robbins problem, IMA J. 
%    Numer. Anal., vol. 13, no. 2, 1993, pp. 215–233.
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
% Date: 28-4-2016, 11:54 (The Netherlands)
% ---------------------------------------------------------
h = 1/n;
r = q*h/2;
r1 = (-1 - r)*ones(1,n-1); r2 = (-1 + r)*ones(1,n-1);
T = diag(2*ones(1,n)) + diag(r1,-1) + diag(r2,1); T = sparse(T);
I = speye(n);
en = I(:,n);
x = h:h:1; y = h:h:1;
a = 2^(alpha -1)*sin(y'); c = 2^(alpha -1)*exp(x');
D1 = spdiags(a,0,n,n);  D2 = spdiags(c,0,n,n);
A = kron(T,I) + kron(I,T) - 2*h*kron(D1,en*en') - 2*h*kron(en*en',D2);
m = n^2;
% In fact, the given exact solution depends on your own choice, it also
% means that the corresponding right hand side vector will be chosen. 
x1 = (-1 + (1-(-1))*rand(m,1)); % generate the random exact solutions.
%x1 = ones(m,1);
b = A*x1;