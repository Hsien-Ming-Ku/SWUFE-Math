function [A,b] = Stecil_Poison(m,opt)
% This function is used to generat the discretized systems from Poisson
% equation in the domain [0,1]
% \times [0,1]\times [0,1], 
%   u_xx + u_yy = 2*(3*x+x^2+y^2),
%   u(x,y) = x^2*(x+y^2)+2 % boundary,
%   u(x,y) = x^2*(x+y^2)+2 % exact solution
% Here we give some comments:
% m+1 means dividing [0,1] into m+1 parts.
% opt means the selection of correction methods.
% by the finite difference method. We can obtain the coefficient matrix A
% and right-hand side vector B. This test problem in modified from [1].
% Reference：
% 1.李厚彪,刘兴平,谷同祥,黄廷祝,李红.Poisson方程有限差分逼近的两种保对称Stencil
%   消元格式[J].计算物理,2010,27(03):335-341.

h = 1 / (m+1);
%Im = speye(m,m);

if opt == 1 %将(i,j+1)处的分量转移到相对格点的对应分量上
   a = [3,11/4 * ones(m-2,1)',3]';
   b = [11/4,2 * ones(m-2,1)',11/4]';
else %将(i,j+1)处的分量转移到临近格点的对应分量上
   a = [3,21/8 * ones(m-2,1)',3]';
   b = [21/8,2 * ones(m-2,1)',21/8]'; 
end
c = ones(m,1);
d = zeros(m,1);d(1) = 1;d(m) = 1;
e = ones(m,1);e(1) = 0;e(m) = 0;
I1 = spdiags([c 0*c c],[-1 0 1],m,m);
I2 = spdiags([c 0*c c],[-2 0 2],m,m);
I3 = spdiags(d,0,m,m);
I4 = spdiags(e,0,m,m);
R = spdiags([-1/4*c a -1/4*c],[-2 0 2],m,m)/h^2;
S = spdiags([-1/8*c 0*c -1/8*c],[-2 0 2],m,m)/h^2;
T = spdiags([-1/8*c -1/4*c -1/8*c],[-1 0 1],m,m)/h^2;
V = spdiags([-1/4*c b -1/4*c],[-2 0 2],m,m)/h^2;
A = kron(I3,R)+kron(I4,V)+kron(I1,S)+kron(I2,T);

%right-hand side vector B
%x = linspace(h,1 - h,m);y=linspace(h,1 - h,m);
%U = sparse(kron(d,d));
%for i = 1:m^2
    %j = ceil(i/m);k = i-(j-1)*m;
    %U(i) = x(j)^2*(x(j)+y(k)^2)+2;
%end % exact solution
% In fact, the given exact solution depends on your own choice, it also
% means that the corresponding right hand side vector will be chosen.
x = rand(size(A,1),1); 
b = A*x; % right-hand side vector b

%spy(A) %画矩阵散点图
