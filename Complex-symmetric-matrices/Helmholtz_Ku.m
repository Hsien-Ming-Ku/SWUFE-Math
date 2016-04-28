%---------------------------------------------------------------
%
%	Helmholtz Equation.
%	
%	Author: Yimin Zhong
%	Last Modified Date: July 5th, 2012
%	Version: 0.0.2
%	Platform: matlab 2009b
%
%---------------------------------------------------------------                                                                      
function [U, H, L1, f1] = Helmholtz_Ku
% Usage:Solve Helmholtz Equation in 2-D square S = [0,1]X[0,1].
% \Delta u + k^2 *u + i*k*\sigma u = 0  in S
% u = f on boundary.   on \partial S
% input: 
%   k: wave number, can be complex number
%   m,n : nodes on the sides.
%set up mesh. 
k = 3;
m = 61;
n = 61;

[x,y] = ndgrid((0:(m-1))/(m-1),(0:(n-1))/(n-1)); 
p = [x(:),y(:)];
t = [1,2,m+2; 1,m+2,m+1]; 
t = kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');
t = kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m); 
b = [1:m,m+1:m:m*n,2*m:m:m*n,m*n-m+2:m*n-1]; 
% b = numbers of all 2m+2n **boundary nodes** preparing for U(b)=0

N = size(p,1);T=size(t,1); % number of nodes, number of triangles
K = sparse(N,N); % zero matrix in sparse format: zeros(N) would be "dense"
F = zeros(N,1); % load vector F to hold integrals of phi's times load f(x,y)
M = sparse(N,N);% mass matrix
U = zeros(N,1);%numerical solution
V = zeros(N,1);
h = min(1/(n-1),1/(m-1)); %mesh size
%sigma
% sigma=0.12*((p(:,1)-0.25).^2+(p(:,2)-0.5).^2<0.01)+0.08*((p(:,1)<0.75).*(p(:,1)>0.5).*(p(:,2)<0.75).*(p(:,2)>0.25))+...
%     0.04*(1-((p(:,1)-0.25).^2+(p(:,2)-0.5).^2<0.01)).*(1-((p(:,1)<0.75).*(p(:,1)>0.5).*(p(:,2)<0.75).*(p(:,2)>0.25)));
sigma = sqrt(3)*ones(N,1);
%construct mass matrix and stiff matrix
for e = 1: T  
  nodes = t(e,:); 
  Pe = [ones(3,1), p(nodes,:)]; 
  Area = abs(det(Pe))/2;
  C = inv(Pe); 
  grad = C(2:3,:);  Ke = Area*(grad'*grad); 
  Me = Area*(1/12)*[[2 1 1];[1 2 1 ];[1 1 2]];
  K(nodes,nodes) = K(nodes,nodes)+Ke; 
  M(nodes,nodes) = M(nodes,nodes)+Me;
end   
%Set up boundary condition
 U(b) = 0.5*exp(1i*(sqrt(5)+1i).*p(b,1,:))+0.5*exp(-1i*(sqrt(5)+1i).*p(b,1,:));
 i = sqrt(-1);
 F = -(K - k^2*M - i*k*diag(sigma)*M)*U;
 freenodes = setdiff(1:N,b); 
 L1 = K(freenodes,freenodes) - k^2*M(freenodes,freenodes) - ...
     i*k*diag(sigma(freenodes))*M(freenodes,freenodes);
 
 f1 = F(freenodes);
 
 [U(freenodes),~] = gmres(L1,f1,30,1e-9);
 H = sigma(:).*U.*conj(U);

 
figure(1)
trisurf(t,p(:,1),p(:,2),0*p(:,1),real(H),'edgecolor','none','facecolor','interp');
view(2),axis equal;colorbar,title('numerical H')
disp('error:');
U_true = 0.5*exp(1i*(sqrt(5)+1i).*p(:,1))+0.5*exp(-1i*(sqrt(5)+1i).*p(:,1));
disp(norm(U-U_true)/norm(U_true));
figure(2)
subplot(2,2,1);
trisurf(t,p(:,1),p(:,2),0*p(:,1),real(U),'edgecolor','none','facecolor','interp');
view(2),axis equal;colorbar,title('numerical real part')
subplot(2,2,2);
trisurf(t,p(:,1),p(:,2),0*p(:,1),imag(U),'edgecolor','none','facecolor','interp');
view(2),axis equal;colorbar,title('numerical image part')
subplot(2,2,3);
trisurf(t,p(:,1),p(:,2),real(U));
subplot(2,2,4);
trisurf(t,p(:,1),p(:,2),imag(U));