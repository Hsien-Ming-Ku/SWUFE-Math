function [A,b] = CCD_Helmholtz1D(n,k,opt)
% This function is used to generate the saddle-point type linear systems
%        [A11   A12 A13]   [1/h^2 *f]   [ 0 ]
%        [A21   A22 A23] * [1/h *f' ] = [ 0 ]
%        [eta*I  0   I ]   [  f''   ]   [ g ],
% Moreover, this system is from the Combined Compact Difference Method for 
% one-dimension Helmholtz equation, 
%          f''(x) + k^2 *f(x) = 1, x\in [0,1]
%                        f(0) = 0,
%            f'(1) - i*k*f(1) = 1;
% refer to [1, pp. 11-32] for details.
% ---------------------------------------------------------------------
% Reference:
% 1. J. Liu, Efficient Preconditioners for the Helmholtz Equation 
%    Discretized by Combined Compact Difference Method, Master thesis, 
%    School of Mathematical Sciences, South China Normal University, May 
%    2010, 70 pages. The file is also available online at 
%    \url{http://lagrange.math.siu.edu/liu/papers/MAThesis_JunLiu.pdf}.
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
% Date: 2-5-2016, 17:42 (The Netherlands)
% --------------------------------------------------------------------
h = 1/(n - 1);
eta = (k*h)^2;
% ---------  Generate the sub-matrix A11  ---------------%
a2 = [0 3*ones(1,n-2)]; a3 = a2(end:-1:1);
if (opt==1 || opt == 3)
   a1 = [-1 -6*ones(1,n-2) -1];
   A11 = diag(a1) + diag(a2,1) + diag(a3,-1);
else
   a1 = [-1 -6*ones(1,n-2) -1i*k];
   A11 = diag(a1) + diag(a2,1) + diag(a3,-1);
end
% ---------  Generate the sub-matrix A12  ---------------%
b1 = [0 (-9/8)*ones(1,n-2)]; b2 = -b1(end:-1:1);
if (opt==1 || opt == 3)
   A12 = diag(b1,1) + diag(b2,-1);
else
   A12 = diag(b1,1) + diag(b2,-1); A12(n,n) = 1/h;
end
% ---------  Generate the sub-matrix A13  -------------- %
c1 = [0 -ones(1,n-2) 0];
c2 = [0 (1/8)*ones(1,n-2)]; c3 = c2(end:-1:1);
A13 = diag(c1) + diag(c2,1) + diag(c3,-1);
% ---------  Generate the sub-matrix A21  --------------- %
d1 = [31 zeros(1,n-2) -31];
d2 = [-32 -15*ones(1,n-2)]; d3 = -d2(end:-1:1);
A21 = diag(d1) + diag(d2,1) + diag(d3,-1); A21(1,3) = 1; A21(n,n-2) = -1;
% ---------  Generate the sub-matrix A22  --------------- %
e1 = [14 16*ones(1,n-2) 14];
e2 = [16 7*ones(1,n-2)]; e3 = e2(end:-1:1);
A22 = diag(e1) + diag(e2,1) + diag(e3,-1);
% ---------  Generate the sub-matrix A23 --------------- %
f1 = [2 zeros(1,n-2) -2];
f2 = [-4 -ones(1,n-2)]; f3 = -f2(end:-1:1);
A23 = diag(f1) + diag(f2,1) + diag(f3,-1);
% ---------- Assemble the global coefficient matrix ----- %
I = speye(n); O = zeros(n,n);
A = [A11 A12 A13;A21 A22 A23;eta*I O I];
A = sparse(A);
% ---------- Assemble the right-hand side b  ------------ %
b = zeros(3*n,1);
xx = (0:h:1).';  % gird nodes;
if opt ==1
   g = (1-pi^2)*(k^2*sin(k*pi*xx));
   b(2*n+1:end) = g;
elseif opt==2
       g = ones(n,1);
       b(2*n+1:end) = g;
elseif opt==3
       g =  (-1+(1-(-1))*rand(n,1));
       b(2*n+1:end) = g;
end
