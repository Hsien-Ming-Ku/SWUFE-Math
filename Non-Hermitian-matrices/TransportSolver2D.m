function TransportSolver2D()
	epsilon=0.1; % diffusion parameter 
	[p,e,t]=initmesh('squareg','hmax',0.05); % mesh 
	np=size(p,2); % number of nodes 
	[A,~,b]=assema(p,t,1,0,1); % diffusion matrix A
					% load vector b 
	bx=ones(np,1); 
	by=ones(np,1); % convection field
	C=ConvectionAssembler2D(p,t,bx,by); % convection matrix C 
       fixed=unique([e(1,:) e(2,:)]); % boundary nodes

	free=setdiff([1:np],fixed); % interior nodes
	b=b(free); % modify b for BC
	A=A(free,free); C=C(free,free); % modify A and C for BC 
	xi=zeros(np,1); % solution vector 
	xi(free)=(epsilon*A+C)\b; % solve for free node values 
	pdesurf(p,t,xi) % plot u
%%%%---------------------------------
function [area,b,c] = HatGradients(x,y)
area=polyarea(x,y);
b=[y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/area;
c=[x(3)-x(2); x(1)-x(3); x(2)-x(1)]/2/area;
%%%%%%%%%%%%%%%%%%%
function C = ConvectionAssembler2D(p,t,bx,by)
np = size(p,2);
nt = size(t,2);
C = sparse(np,np);
for i=1:nt
    loc2glb = t(1:3,i);
    x = p(1,loc2glb);
    y = p(2,loc2glb);
    [area,b,c] = HatGradients(x,y);
    bxmid = mean(bx(loc2glb));
    bymid = mean(by(loc2glb));
    CK = ones(3,1)*(bxmid*b+bymid*c)'*area/3;
    C(loc2glb,loc2glb)=C(loc2glb,loc2glb)+CK;
end