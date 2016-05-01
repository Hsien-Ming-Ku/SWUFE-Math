function [A,b] = Convection_diffusionFEM2D
% solves the 2D convection-diffusion equation
%          -nu (u_xx+u_yy) + wx u_x + wy u_y + la u = f
% in the square domain [0,1] x [0,1] with Dirichlet boundary conditions 
% u = g using piecewise linear triangular finite elements; refer to 
% [1, pp. 130-134] and [2, pp. 55-63, 76-79.]
% ------------------------------------------------------------
% Reference
% 1. A. Greenbaum, Iterative Methods for Solving Linear Systems, SIAM, 
%    Philadelphia, USA, 1997.
% 2. Y. Saad, Iterative Methods for Sparse Linear Systems, 2nd edition, 
%    SIAM, Philadelphia, USA, 2003.
% ---------------------------------------------------------
% Modified by: 
% Mr. Xian-Ming Gu, who works as a PhD student in School of Mathematical
% Sciences, University of Electronic Science and Technology of China 
% (UESTC). From 1 August, 2014 to 1 Augest, 2016, He also work in the 
% Johann Bernoulli Institute for Mathematics and Computer Science (JBI), 
% University of Groningen (RUG). Prof. Ting-Zhu Huang (UESTC) and Dr. Bruno
% Carpentieri (RUG) are Xian-Ming's supervisor for PhD program.
% --------------------------------------------------------------------
% Email: guxianming@live.cn
% URL: https://github.com/Hsien-Ming-Ku
% Date: 28-4-2016, 13:30 (The Netherlands)
% ---------------------------------------------------------
% define pb
param = [1  ...       % nu (viscosity)
         10  ...       % wx (x-convection)
         -20  ...       % wy (y-convection)
         1];          % lambda (linear term)
%-----------------------------------------------------------------------
n = 80;               % pb size n^2
x = [0;rand(n,1);1];  % define random finite element mesh
y = [0;rand(n,1);1];
x = [zeros(n+2,1);x;x;ones(n+2,1);rand(n^2,1)];
y = [y;zeros(n+2,1);ones(n+2,1);y;rand(n^2,1)];
bdy = [ones(n+2,1);ones(n+2,1);ones(n+2,1);ones(n+2,1);Inf*ones(n^2,1)];
T = delaunay(x,y);    % T is nt x 3. Each row defines 3 vertices.
nv = length(x);
nt = size(T,1);
figure(2);
subplot(1,2,1)
trimesh(T,x,y,zeros(nv,1),'EdgeColor',[0 0 0]) % show triangulation
view(2)
axis square
axis off
title('Finite element mesh')
%-----------------------------------------------------------------------
[A,b] = system(T,x,y,bdy,param,@f,@g); % build linear system
figure(1)
subplot(1,2,1)
permute = 1;
switch permute
case 0 % no permutation
    p = 1:size(A,1); type = '';
case 1 % minimum degree algorithm
    p = symamd(A); type = ' (Min Degree ordering)';
case 2 % reverse Cuthill-McKee
    p = symrcm(A); ' (RCM ordering)';
end
spy(A(p,p))
axis square
set(gca,'XTick',[],'YTick',[])
title(['Sparsity pattern' type])
nz = nnz(A);
xlabel(['NNZ = ' int2str(nz) '    density = ' num2str(nz/nv^2)])
subplot(1,2,2)
[L,U] = lu(A(p,p));
spy(L+U)
title('Sparsity pattern of factorization (L+U)')
nz = nnz(L+U);
axis square
set(gca,'XTick',[],'YTick',[])
xlabel(['NNZ = ' int2str(nz) '    density = ' num2str(nz/nv^2)])
%-----------------------------------------------------------------------
i = find(bdy==Inf); 
j = find(bdy==1);
u(i(p)) = A(p,p)\b(p);      % solve for interior values
u(j) = feval(@g,x(j),y(j)); % boundary values
figure(2)
subplot(1,2,2)
ts = trisurf(T,x,y,u);
axis square
shading interp
set(ts,'EdgeColor',[0 0 0])
xlabel('x','FontWeight','bold')
ylabel('y','FontWeight','bold')
zlabel('u(x,y)','FontWeight','bold')
set(gca,'XTick',[0 1],'YTick',[0 1])
title('Numerical solution')
%=======================================================================
function [A,b] = system(T,x,y,bdy,param,f,g)
% generate finite element system
npts = 3;
switch npts
case 3 % 3 points quadrature rule
    c = [4 1 1;1 4 1;1 1 4]/6; 
    wt = [1 1 1]'/3;
case 6 % 6 points quadrature rule
    alfa = 0.091576213509771;
    beta = 0.445948490915965;
    c = [1-2*alfa alfa alfa;alfa 1-2*alfa alfa;alfa alfa 1-2*alfa; ...
            1-2*beta beta beta;beta 1-2*beta beta;beta beta 1-2*beta];
    w1 = 0.109951743655332;
    w2 = 0.223381589678011;
    wt = [w1 w1 w1 w2 w2 w2]';
otherwise
    disp('quadrature rule not implemented')
end
I3 = ones(3,1); 
M = [2 1 1;1 2 1;1 1 2]/12;
nt = size(T,1);
nv = length(x);
xt = reshape(x(T),nt,3)';
yt = reshape(y(T),nt,3)';
dx = diff(xt); 
dy = diff(yt); 
d = dx(1,:).*dy(2,:)-dx(2,:).*dy(1,:);
d = kron(d,I3);
area = abs(d)/2;
gx = (yt([2 3 1],:)-yt([3 1 2],:))./d; % basis function x-gradient
gy = (xt([3 1 2],:)-xt([2 3 1],:))./d; %                y-gradient
cx = kron(gx,I3); 
cy = kron(gy,I3);
lapl = cx.*kron(I3,gx)+cy.*kron(I3,gy);
M = kron(M(:),ones(1,nt));
all = param(1)*lapl+(param(2)/3)*cx+(param(3)/3)*cy+param(4)*M;
all = all.*kron(area,I3);
b = area.*(c'*(kron(wt,ones(1,nt)).*feval(f,c*xt,c*yt)));
T = T';
i = kron(I3,T);
i = i(:);
j = kron(T(:),I3);
A = sparse(i,j,all(:));
i = T(:);
j = ones(3*nt,1);
b = sparse(i,j,b(:));
b = full(b);
i = find(bdy==Inf); % Dirichlet boundary conditions                  
j = find(bdy==1);                     
ub = feval(g,x(j),y(j));                          
b = b(i)-A(i,j)*ub; 
A = A(i,i); 
%=======================================================================
function f = f(x,y)
f = 100*(1-x).^2.*y;
%=======================================================================
function g = g(x,y)
g = x.^2.*y;
