function [A,b] = cavity_Helmholtz(q,omega,theta)
% Let us consider the electromagnetic scattering problem from a large 
% rectangular cavity on the x-y plane in which the medium is y-directional 
% inhomogeneous. In the transverse magnetic polarization case, when the 
% model Helmholtz equation with positive wave number is discretized by 
% the five-point finite difference scheme with uniform stepsize h, we 
% obtain a block two-by-two system of linear equations with following form
%   A*x = b: 
%               A = [ B  E ]
%                   [ F  C ]
% q+1 means dividing [0,1] into q+1 parts.
% omega is a real constant
% theta ≥ 0 is a real constant
% Reference
% 1. Z.-Z. Bai, Structured preconditioners for nonsingular matrices of 
%    block two-by-two structures, Math. Comp., 75 (2006), pp. 791-815.
% ------------------------------------------------------------------
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
% Date: 28-4-2016, 12:06 (The Netherlands)
% ---------------------------------------------------------^
h = 1/(q+1);
% p = q^2; 
I = speye(q);
eq = [zeros(q-1,1); 1];
r = ones(q,1);
V = spdiags([(-1+1/2*theta*h)*r 2*r (-1-1/2*theta*h)*r],[-1 0 1],q,q);
Omega = spdiags((h^2)*omega.^2,0,q,q);
B = kron(V,I) + kron(I,V) - kron(I,Omega);
% B is symmetric positive definite
G = zeros(q,q);
for i = 1:q
    for j = 1:q
        G(i,j) = 1/(i+j)^2;
    end
end
%I1 = speye(p);
C = I - h*G; 
% C is symmetric positive/negative definite.
E = kron(I,eq); F = -E';
A = [B E;F C];
% In fact, the given exact solution depends on your own choice, it also
% means that the corresponding right hand side vector will be chosen. 
% -------------------------------------------------------------------
x_ext =(-1+(1-(-1))*rand(size(A,1),1))+1i*(-1+(1-(-1))*rand(size(A,1),1);
%x_ext = (-1+(1-(-1))*rand(size(A,1),1));
%x_ext = ones(size(A,1),1);
b = A*x_ext; 
