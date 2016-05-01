% ExPDE21
% What does convection look like?  Reduced diffusion
% classic method

% u_t = Du_xx + au_x + cu + S(x,t)


dx = 1e-2;
dt = 1e-5;
x = 0:dx:1;
imax = length(x);
Tmax = .5;
nmax = round(Tmax/dt) + 1;

D = .01;
b = -3;  % + bu_x
a = 0;  % + au
S = @(x,t) 0;
f = @(t) 0; % left BC
g = @(t) 0; % right BC
h = @(x) exp(-100*(x-.4).^2); % IC

rt = D*dt/dx^2;
beta = b*dt/(2*dx);
np = 5*(nmax-1)/100;

u = zeros(imax,2);

% initial condition
u(:,1) = h(x);
figure(1)
plot(x,u(:,1),'r')
hold on

for n=1:nmax-1
    tn = (n-1)*dt;
    u(1,2) = f(tn+dt);
    for i=2:imax-1
        uxx = rt*(u(i-1,1)-2*u(i,1)+u(i+1,1));
        ux = beta*(u(i+1,1)-u(i-1,1));
        u(i,2) = u(i,1) + uxx + ux + dt*a*u(i,1) + dt*S(x(i),tn);
    end
    u(imax,2) = g(tn+dt);
    if (mod(n,np)==0)
        plot(x,u(:,2),'b')
        hold on
        pause
    end
    u(:,1) = u(:,2);
end
