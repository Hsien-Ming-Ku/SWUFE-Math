function [A, b, x_ext] = Helmholtz_Robbins(m, k, opt)
% This function is used to generate the nonsymmetric linear systems A*x = b
% coming from the finite difference discretization of 2D Helmhotz problem 
% with Robbins boundary conditions in square domian [0,1]^2:
%   -(u_xx + u_yy) - k^2 * u = f,
%      u(0, y) = g(y), u(x, 0) = \rho(x);
%      u_x(1, y) =  p*u(1, y) + a(y);
%      u_y(x, 1) = q*u(x, 1) + c(x);
% where k, p and q are constants with finite difference grid size h = 1/m; 
% refer to [1] for details.
% here:
% 1) m: is the number of grid nodes in x-direction (or y-direction);
% 2) k: is the wave-number;
% 3) opt: four kinds of different test problems.
% Output:
% A:     the coefficient matrix of discretized linear systems;
% b:     the right-hand side vector of discretized linear systems;
% x_ext: the discretized vector of original exact solutions.
% ---------------------------------------------------------------------
% Reference
% 1. L. Jiang, J. Huang , X.-G. Lv, Q.-S. Cheng, A preconditioned method 
%    for the solution of the Robbins problem for the Helmholtz equation, 
%    ANZIAM J., 52 (2010), pp. 87-100.
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
% Date: 23-5-2016, 14:41 (The Netherlands)
% ----------------------- Physics parameters ---------------------
h = 1/m;  % grid size
if opt==1
   p = 1; q = 1/2;  % opt = 1; Example 4.1
elseif opt==2
   p = -1; q = 1;   % opt = 2; Example 4.2
elseif opt==3
   p = 1; q = -1;   % opt = 3; Example 4.3
elseif opt==4
   p = -1; q = -1;  % opt = 4; Example 4.4 (not show in Ref. [1])
%elseif opt==5
%  p = ?; q = ?   % depends on your decision
end
lamda = k^2;  % lamda = k^2, k is the wave-number
Im = speye(m);  % Identity matrix with order m
% ----- Generate the coefficient (tridiagonal) matrix B in Ref. [1] ------
e1 = [(4 - lamda*h^2)*ones(m - 1,1); 4 - lamda*h^2 - 2*h*q];
e2 = [ones(m-2,1);2;1];
B = spdiags([-e2 e1 -ones(m,1)],[-1 0 1],m,m);
% ----- Generate the coefficient (tridiagonal) matrix C in Ref. [1] ------
ex = [zeros(m-1,1);-2*h*p];
C = spdiags([-e2 ex -ones(m,1)],[-1 0 1],m,m);
% --- Generate the coefficient (bloc tridiagonal) matrix A in Ref. [1] ---
A = kron(Im,B) + kron(C,Im);
% Step 2 assemble the right hand side vector b and exact solution u.
if opt==1   % we set: u(x,y) = e^(xy)   f(x,y) = -(x^2 + y^2 + k^2)e^(xy)
   for k = 1:m
       for j = 1:m
           u(k,j) = exp(k*h*j*h);    % compute the exact solution vector
           f(k,j) = -(((k*h)^2 + (j*h)^2 + ...
                    lamda)*u(k,j))*h^2; % Corresponding to the right hand 
                                      % side vector from f at inner point.
       end 
   end
   for k = 1:m
       f(1,k) = f(1,k) + 1;   f(k,1) = f(k,1) + 1;
       c(k)= (k*h - q)*exp(k*h);  f(k,m) = f(k,m) + 2*h*c(k);
       a(k) = (k*h - p)*exp(k*h);  f(m,k) = f(m,k) + 2*h*a(k);
   end
   f = f';
   for k = 1:m
       b(m*(k-1)+1:m*k,1:1) = f(1:m,k:k);
   end
   u = u';
   for k = 1:m
       uvector(m*(k-1)+1:m*k,1:1) = u(1:m,k:k);
   end
elseif opt==2 %  u(x,y) = sin(Pi*x/2)sin(Pi*y),f(x,y)=(-5Pi^2/4-lamda)u
    for k = 1:m
        for j = 1:m
            u(k,j) = sin(pi*k*h/2)*sin(pi*j*h);
            f(k,j) = (5*pi*pi/4 - lamda)*u(k,j)*h^2;
        end 
    end
    for k=1:m
        c(k) = -pi*sin(pi*k*h/2);   f(k,m) = f(k,m) + 2*h*c(k);
        a(k) = sin(pi*k*h);   f(m,k) = f(m,k) + 2*h*a(k);
    end
    f=f';
    for k=1:m
        b(m*(k-1)+1:m*k,1:1) = f(1:m,k:k);
    end
    u = u';
    for k=1:m
        uvector(m*(k-1)+1:m*k,1:1)=u(1:m,k:k);
    end
elseif opt==3  % u(x,y) = x^2 + y^2   f(x,y) = -(4 + lamda*(x^2+y^2))
    for k = 1:m
        for j = 1:m
            u(k,j) = (k*h)^2 + (j*h)^2;
            f(k,j) = (-4 - lamda*u(k,j))*h^2;
        end 
    end
    for k = 1:m
        f(1,k) = f(1,k) + (k*h)^2;  f(k,1) = f(k,1) + (k*h)^2;
        c(k) = 3 + (k*h)^2;   f(k,m) = f(k,m) + 2*h*c(k);
        a(k) = 1 - (k*h)^2;   f(m,k) = f(m,k) + 2*h*a(k);
    end
    f = f';
    for k = 1:m
        b(m*(k-1)+1:m*k,1:1) = f(1:m,k:k);
    end
    u = u';
    for k = 1:m
        uvector(m*(k-1)+1:m*k,1:1) = u(1:m,k:k);
    end
elseif opt==4   % u(x,y) = e^x + e^y  
    for k=1:m
        for j=1:m
            u(k,j)=exp(k*h+j*h);
            f(k,j)=(-2 - lamda)*u(k,j)*h^2;
        end 
    end
    for k=1:m
        f(1,k) = f(1,k) + exp(k*h);   f(k,1) = f(k,1) + exp(k*h);
    end
    for k=1:m
        c(k) = 2*exp(k*h + 1);  f(k,m) = f(k,m) + 2*h*c(k);
        a(k) = 2*exp(k*h+1);    f(m,k)=f(m,k)+2*h*a(k);
    end
    f = f';
    for k=1:m
        b(m*(k-1)+1:m*k,1:1) = f(1:m,k:k);
    end
    u = u';
    for k=1:m
        uvector(m*(k-1)+1:m*k,1:1)=u(1:m,k:k);
    end
end
x2 = A\b;
x_ext = uvector;
max_res = norm(uvector - x2,inf)/norm(uvector,inf);
%disp(max_res);
fprintf('-- L_{inf}-norm error: %.10e\n ', max_res); % Error.