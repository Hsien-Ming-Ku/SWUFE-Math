function [L,U] = indefiniteILU(A)
% This function is to handle the ILU decomposition for (indefinite) matrix,
% it means that there are some zeros in its diagonal part, this situation
% is not suitable for using ILU decomposition directly. We give a simple
% strategy [1-2] for remedying this difficult. We denote the matrix as
%              A = (a_ij) \in mathcal{C}^{n\times n}
% then we obtain the A1 = A + epst*max(|a_ii|)*I_n, then we compute ILU
% decomposition of A1 for providing the ILU preconditioner for original
% coefficient matrix A.
% ---------------------------------------------------------------------
% References:
% 1. E. Chow and Y. Saad, Experimental study of ILU preconditioners for 
%    indefinite matrices, J. Comput. Appl. Math., 86 (1997), pp. 387-414.
% 2. S. Güttel and J. Pestana, Some observations on weighted GMRES, Numer. 
%    Algorithms, 67 (2014), pp. 733-752.
% ---------------------------------------------------------------------
% Written by: 
% Mr. Xian-Ming Gu, who works as a PhD student in School of Mathematical
% Sciences, University of Electronic Science and Technology of China 
% (UESTC). From 1 August, 2014 to 1 Augest, 2016, He also work in the 
% Johann Bernoulli Institute for Mathematics and Computer Science (JBI), 
% University of Groningen (RUG). Prof. Ting-Zhu Huang (UESTC) and Dr. Bruno
% Carpentieri (RUG) are Xian-Ming's supervisor for PhD program.
% --------------------------------------------------------------------
% Email: guxianming@live.cn
% URL: https://github.com/Hsien-Ming-Ku
% Date: 4-6-2016, 00:43 (The Netherlands)
% --------------------------------------------------------
d = diag(A);
n = size(A,1);
t = max(abs(d));
sr = min(abs(d));
if t==0   % all diagonal entries are zeros
   sigma = 1e-12; % recommended, even you can try "simga = 1e-14"
elseif (t~=0 && sr==0) % some but not all diagonal elements a_ii of A zero
   sigma = (1e-12)*t;  % recommended, even you can try "simga = (1e-14)*t"
else
   sigma = 0; % no zeros in its diagonal part.  
end
% ------------------ creat the auxiliary matrix A1 ----------------------
A1 = A + sigma*speye(n,n);
% ------- choose the suitable kind of ILU decomposition for A1 -----------
setup.type = 'crout';  % can be revised for your study, 
setup.milu = 'row';    % can be revised for your study, 
setup.droptol = 0.01;  % can be revised for your study, e.g., 0.001;
% the above choice is just an example !! you can obtain the suitable ILU
% decomposition by your own decision. Please refer to "ilu doc" !!!
% ------------------------------------
[L,U] = ilu(A1,setup); % Then L, U is just the targeted ILU preconditioner
                       % for original coefficient matrix A.
clear A1;



