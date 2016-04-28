function [A,b] = convection_diffusionHOC2D(n,tau,opt)
% This function is to generate the linear systems A*x = b arising from the
% fourth-order difference scheme for 2D (steady) convection-diffusion 
% equation, refer to [1] for details. 
%%%%%%%%%%%%%
% References
% 1. M.S. Sunhaloo, R. Boojhawon, A. Gopaul, M. Bhuruth,On block-circulant 
%    preconditioners for high-order compact approximations of convection
%    diffusion problems, J. Comput. Appl. Math., 234 (2010), pp. 1312-1323.
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
% Date: 28-4-2016, 11:41 (The Netherlands)
% -------------------------------------------------------
h = 1/(n + 1);
gamma = tau*h/2;
r1 = [-(4 + 4*gamma + 2*gamma^2) -(1 + gamma) zeros(1,n - 2)];
L1 = sparse(toeplitz(r1));
r2 = [20+ 4*gamma^2 -4 zeros(1,n-2)];
K1 = sparse(toeplitz(r2));
r3 = [-(4 - 4*gamma + 2*gamma^2) -(1 - gamma) zeros(1,n - 2)];
M1 = sparse(toeplitz(r3));
I1 = sparse(diag(ones(1,n))); I2 = sparse(diag(ones(1,n - 1),-1));
I3 = sparse(diag(ones(1,n - 1),1));
A = kron(I1,K1) + kron(I2,L1) + kron(I3,M1); % Block tri-diagonal matrix.
if opt==1
   b = A*ones(length(A),1);
else
   b = A*rand(length(A),1);
% Here the right-hand side vector can be chosen by yourself.   
end


