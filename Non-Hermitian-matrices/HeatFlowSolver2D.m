function [Afem,r] = HeatFlowSolver2D()
% Here we employ the function "HeatFlowSolver2D" generate a linear system 
% corresponding to the stabilized convection-diffusion or transport 
% equation (2D), here in the weak form:
%           a(u, v)  =  \ell(v);
% a(u,v) = \nu(grad u, grad v) +(w*grad u, v) + \delta(w*grad u, w*grad v)
% \ell(v) = (f, v);
% refer to [1].
% ------------------------------------------------------------
% Reference
% M. G. Larson and F. Bengzon, Transport Problems, in The Finite Element 
% Method: Theory, Implementation, and Applications, volume 10 of Texts in 
% Computational Science and Engineering, Springer, 2013, pp. 241-256.
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
% Date: 2-5-2016, 00:43 (The Netherlands)
% ---------------------------------------------------------
	%channel=RectCircg(); % channel geometry
    channel=[2	2	2	2 	1	1	1	1
	6	6	-2	-2	-1	0	1	0
	6	-2	-2	6	0	1	0	-1
	-2	2	-2	-2	0	-1	0	1
	2	2	2	-2	-1	0	1	0
	1	1	0	1	0	0	0	0
	0	0	1	0	1	1	1	1
	0	0	0	0	0	0	0	0 
	0	0	0	0	0	0	0	0 
	0	0	0	0	1	1	1	1];
	epsilon = 0.01; % diffusion parameter
	h = 0.05; % mesh size
	[p,e,t]=initmesh(channel,'hmax',h); % create mesh
    pdemesh(p,e,t)
    pause
    % pdemesh(p,e,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	A=assema(p,t,1,0,0); % stiffness matrix
	x=p(1,:); y=p(2,:); % node coordinates 
    [bx,by]=FlowField(x,y); % evaluate vector field b 
 	C=ConvectionAssembler2D(p,t,bx,by); % convection matrix 
 	Sd=SDAssembler2D(p,t,bx,by); % GLS stabilization matrix 
	[R,r]=RobinAssembler2D(p,e,@Kappa2,@gD2,@gN2); % Robin BC 
    delta=h; % stabilization parameter
    Afem = epsilon*A-C'+R+delta*Sd;
	U = Afem\r; % solve linear system 
	pdecont(p,t,U), axis equal % plot solution
    colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g = RectCircg()
	g = [2	2	2	2 	1	1	1	1
	6	6	-2	-2	-1	0	1	0
	6	-2	-2	6	0	1	0	-1
	-2	2	-2	-2	0	-1	0	1
	2	2	2	-2	-1	0	1	0
	1	1	0	1	0	0	0	0
	0	0	1	0	1	1	1	1
	0	0	0	0	0	0	0	0 
	0	0	0	0	0	0	0	0 
	0	0	0	0	1	1	1	1]; 


function [bx,by] = FlowField(x,y)
    a=1; % cylinder radius
    Uinf=1; % free stream velocity
    r2=x.^2+y.^2; % radius vector squared
    bx=Uinf*(1-a^2*(x.^2-y.^2)./r2.^2); % x-component of b 
    by=-2*a^2*Uinf*x.*y./r2.^2; % y-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = ConvectionAssembler2D(p,t,bx,by) 
	np=size(p,2);
	nt=size(t,2);
	C=sparse(np,np);
	for i=1:nt
		loc2glb=t(1:3,i);
		x=p(1,loc2glb);
		y=p(2,loc2glb); 
		[area,b,c]=HatGradients(x,y); 
		bxmid=mean(bx(loc2glb)); 
		bymid=mean(by(loc2glb)); 
		CK=ones(3,1)*(bxmid*b+bymid*c)'*area/3; 
		C(loc2glb,loc2glb)=C(loc2glb,loc2glb)+CK;
    end
 
    
function [area,b,c] = HatGradients(x,y) 
    area=polyarea(x,y);
    b=[y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/area; 
    c=[x(3)-x(2); x(1)-x(3); x(2)-x(1)]/2/area;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Sd = SDAssembler2D(p,t,bx,by) 
	np=size(p,2);
	nt=size(t,2);
	Sd=sparse(np,np);
    for i=1:nt
        loc2glb=t(1:3,i);
        x=p(1,loc2glb);
        y=p(2,loc2glb);
        [area,b,c]=HatGradients(x,y); 
        bxmid=mean(bx(loc2glb)); 
        bymid=mean(by(loc2glb)); 
        SdK=(bxmid*b+bymid*c)*(bxmid*b+bymid*c)'*area; 
        Sd(loc2glb,loc2glb)=Sd(loc2glb,loc2glb)+SdK;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R,r] = RobinAssembler2D(p,e,kappa,gD,gN) 
	R = RobinMassMatrix2D(p,e,kappa);
	r = RobinLoadVector2D(p,e,kappa,gD,gN);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = RobinMassMatrix2D(p,e,kappa)
	np = size(p,2); % number of nodes
	ne = size(e,2); % number of boundary edges
	R = sparse(np,np); % allocate boundary matrix 
    for E = 1:ne
		loc2glb = e(1:2,E); % boundary nodes
		x = p(1,loc2glb); % node x-coordinates
		y = p(2,loc2glb); % node y-
		len = sqrt((x(1)-x(2))^2+(y(1)-y(2))^2); % edge length 
		xc = mean(x); yc = mean(y); % edge mid-point
		k = kappa(xc,yc); % value of kappa at mid-point
		RE = k/6*[2 1; 1 2]*len; % edge boundary matrix 
		R(loc2glb,loc2glb) = R(loc2glb,loc2glb) + RE;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = RobinLoadVector2D(p,e,kappa,gD,gN) 
    np = size(p,2);
    ne = size(e,2);
    r = zeros(np,1);
    for E = 1:ne
        loc2glb = e(1:2,E);
        x = p(1,loc2glb);
        y = p(2,loc2glb);
        len = sqrt((x(1)-x(2))^2+(y(1)-y(2))^2);
        xc = mean(x); yc = mean(y);
        tmp = kappa(xc,yc)*gD(xc,yc)+gN(xc,yc); 
        rE = tmp*[1; 1]*len/2;
        r(loc2glb) = r(loc2glb) + rE;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
function k = Kappa2(x,y) 
    k=0;
    if x<-1.99 % inflow
    	k=1e6; 
    end
    if sqrt(x^2+y^2)<1.01 % cylinder 
    	k=1e6;
    end
    if x>5.99 % outflow
        [bx,by]=FlowField(x,y);
        nx=1; ny=0; % normal components 
      k=bx*nx+by*ny; % kappa = dot(b,n)
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g = gD2(x,y)
	g=0;
	if sqrt(x^2+y^2)<1.01, 
		g=1; 
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g = gN2(x,y) 
	g=0;

