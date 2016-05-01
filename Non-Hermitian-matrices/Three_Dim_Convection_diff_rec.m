function [A,b] = Three_Dim_Convection_diff_rec(n,q,p,opt)
% This function is to generate the linear systems A*x = b arising from the
% fourth-order difference scheme for 3D (steady) convection-diffusion 
% reaction equation, 
%       -(u_xx + u_yy+ u_zz) + q(u_x + u_y + u_z) - p*u = f,
% refer to [1] for details. 
% Here we note that you can choose larger p for generating ill-conditioned
% matrix.
%%%%%%%%%%%%%
% References
% 1. E. Jiang, Algorithm for solving shifted skew-symmetric linear system,
%    Front. Math. China., 2 (2007), pp. 227-242.
% ----------------------------------------------------------------------
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
% Date: 2-5-2016, 00:25 (The Netherlands)
% -------------------------------------------------------------------
h = 1/(n + 1);
r = q*h/2;
e = ones(n,1);
It = speye(n);
Tx = spdiags([(-1 - r)*e, 6*e, (-1 + r)*e],[-1,0,1],n,n);
Ty = spdiags([(-1 - r)*e, 0*e, (-1 + r)*e],[-1,0,1],n,n);
Tz = spdiags([(-1 - r)*e, 0*e, (-1 + r)*e],[-1,0,1],n,n);

A1 = kron(kron(Tx,It),It) + kron(kron(It,Ty),It) + kron(kron(It,It),Tz);
A = A1 - p*speye(n^3);
if opt==1
   b = A*ones(n^3,1);
else
   b = A*randn(n^3,1);
% Here the right-hand side vector can be chosen by yourself.   
end
