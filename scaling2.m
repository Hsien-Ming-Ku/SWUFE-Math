function [DA,Db,D1] = scaling2(A, b, p, opt)
% The function [DA,Db] = scaling(A,b) implements the
% scaling techniques (including Row scaling and Column scaling) for the
% original linear systems. It often can reduce the condition number of the
% original coefficient matrix. In fact, it can be regarded as a simply
% preconditioning techniques.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
%       A    Original coefficient matrix
%       b    Right hand side vector
%       p    Parameter for the L_p norm, p >= 1, refer to [2]
%       opt  Choices for each kind of scaling techniques, such as Row
%            scaling, two-side scaling;
% Output:
%
%       DA   New matrix after Row (or Column) scaling.
%       Db   Corresponding right hand side vector.
% --------------------------------------------------------------------
% References
% 1. B. Carpentieri, Y.-F. Jing, T.-Z. Huang, The BiCOR and CORS iterative 
%    algorithms for solving nonsymmetric linear systems, SIAM J. Sci. 
%    Comput., 33 (2011), pp. 3020--3036.
% 2. D. Gordon, R. Gordon, Row scaling as a preconditioner for some 
%    nonsymmetric linear systems with discontinuous coefficients, J. 
%    Comput. Appl. Math., 234 (2010), pp. 3480--3495.
% 3. G. Gambolati, G. Pini, M. Ferronato, Scaling improves stability of 
%    preconditioned CG-like solvers for FE consolidation equations, Int. J. 
%    Numer. Anal. Methods Geomech. 27 (2003), pp. 1043--1056.
% 4. A. van der Sluis, Condition numbers and equilibration of matrices, 
%    Numer. Math., 14 (1969), pp 14--23.
% -----------------------------------------------------------------------
% Developped (or copyright) by Hsien-Ming Ku, who is an associate professor
% of the School of Mathematics, Southwestern University of Finance and Economics, Chengdu, 
% 611130, P. R. China. 
% Contact us by e-mail,
% E-mail: guxianming@live.cn
% Thank to Dr. Jing Meng and Prof. Liang Li for their kind suggestions to
% modify the codes.
% Date: 2013-12-16; (Chengdu)
% Posted date: 2-5-2016, 17:55 (the Netherlands)
% updated date: 30/7/2024 (Chengdu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin <2 || isempty(b)
    b = [];
end
% Set up the Scaling process,
A = sparse(A);
n = size(A,1);       % the length of the coefficient matrix.
if opt==1
   r1 = (1./max(abs(A))).'; % the maximum number in each column
   r2 = (1./max(abs(A.'))).'; % the maximum number in each row
   D1 = spdiags(sqrt(r1),0,n,n);
   D2 = spdiags(sqrt(r2),0,n,n);
   DA = D2*(A*D1);
   Db = D2*b;
elseif opt==2
       r3 = (sum((abs(A.')).^p)).^(1/p);
       r3 = r3.';
       D = spdiags(1./r3,0,n,n);
       DA = D*A;
       Db = D*b;
       D1 = [];
elseif opt==3
       r4 = max(abs(A.'));             % the maximum number in each row
       r4 = r4.';
       D = spdiags(1./r4,0,n,n);
       DA = D*A;
       Db = b;  
elseif opt==4
       r1 = (sum((abs(A)).^p)).^(1/p); % the maximum number in each column
       r1 = (1./r1).';
       r2 = (sum((abs(A.')).^p)).^(1/p); % the maximum number in each row
       r2 = (1./r2).';
       D1 = spdiags(sqrt(r1),0,n,n);
       D2 = spdiags(sqrt(r2),0,n,n);
       DA = D2*A*D1;
       Db = D2*b;
end
% End this SCALING2 codes.
