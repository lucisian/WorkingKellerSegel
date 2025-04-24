clear
%set the number of dimensions

dims=1;
%number of gridpoints per dimension

m=4000;
%Total number of gridpoints; varies by dimension as:
%2D make N=m^2; 1D make N=m;
if(dims==1)
N=m;
elseif(dims==2)
N=m^2;
end

%parameters in the reaction kinematics
a=144;
D=1;
p0=(a+D+2*sqrt(a*D));
p2=1;
p4=1;
epsilon=0.1;
p=p0+epsilon*epsilon*p2+epsilon*epsilon*epsilon*epsilon*p4;

%Domain length
L=3*2*pi*sqrt(sqrt(D/a));

%Spatial step size
dx=L/(m-1);

%Time interval to solve the equations on
T=linspace(0,10000,100);

%Spatial domain (needed for plotting only)
x=linspace(0,L,m);

% (Sparse) Laplacian matrix
e=ones(m,1);
Lap=spdiags([e,-2*e,e],[1,0,-1],m,m);
Adv=spdiags([e,-e],[1,-1],m,m);

%Neumann
Lap(1,1)=-1;
Lap(end,end)=-1;
Adv(1,1)=-1;
Adv(end,end)=1;

%Indices corresponding to u variable and v variable. THis lets us stack them both in a vector U and write u = U(ui) and v = U(vi).
ui=1:N;
vi=N+1:2*N;

%Reaction kinetics
f=@(u,v)u.*(1-u);
g=@(u,v)u-a*v;

if(dims==1)
    %1D Laplacian
    Lap = (1/dx)^2*Lap;
    Adv=(1/(2*dx))*Adv;
    %chi = @(U)U./(1+U.^(2));
    F = @(t,U)[f(U(ui),U(vi)) + Lap*U(ui) -...
    p*(U(ui).*(Lap*U(vi)) + (Adv*(U(ui))).*(Adv*U(vi)));
    g(U(ui),U(vi)) + D*Lap*U(vi)];
elseif(dims==2)
    %2d Laplacian
    I = speye(m);
    Advx = (1/(2*dx))*kron(Adv,I);
    Advy = (1/(2*dx))*kron(I, Adv);
    Lap = (1/dx)^2*(kron(Lap,I) + kron(I, Lap));

    %Put together the reaction kinetics+diffusion terms into a big vector
    F = @(t,U)[f(U(ui),U(vi)) + Lap*U(ui) -...
    p*(U(ui).*(Lap*U(vi)) +...
    (Advx*U(ui)).*(Advx*U(vi)) + (Advy*U(ui)).*(Advy*U(vi)));
    g(U(ui),U(vi)) + D*Lap*U(vi)];
end

rng(1);

%Initial condition - this is a small normally distributed perturbation of the homogeneous steady state of our kinetics
U0 = [1*ones(N,1);(1/a)*ones(N,1)]+1e-2*randn(2*N,1);
JacSparse=sparse([Lap,Lap;eye(N),Lap]);
odeOptions=odeset('JPattern',JacSparse);%,'RelTol',1e-5,'AbsTol',1e-5,'maxstep',0.1);
[T,U]=ode15s(F,T,U0,odeOptions);
P=U(end,ui);
close all;

%Third order coefficient - supercritical when positive, subcritical when
%engative
C3 = (9*a-67*sqrt(a*D)-16*D)/72;

if(dims==1)
    plot(x,U(end,ui),'linewidth',2); hold on

    set(gca,'fontsize',24);

    hold on
    
    k=sqrt(sqrt(a/D));
    A0=-1.01751;
    wu=1-A0^(2)*epsilon^(2)/2 +2027*A0^(2)*epsilon^(4)/28730 -23297*A0^(4)*epsilon^(4)/12888720 ...
        +(A0*epsilon +A0*epsilon^(3)/2210 +37*A0^(3)*epsilon^(3)/236720)*cos(k*x)...
        +(50*A0^(2)*epsilon^(2)/27 +954883*A0^(2)*epsilon^(4)/41888340 -20377519009*A0^(4)*epsilon^(4)/11135854080)*cos(2*k*x)...
        +(140525*A0^(3)*epsilon^(3)/55296)*cos(3*k*x)...
        +(137709271*A0^(4)*epsilon^(4)/40310784)*cos(4*k*x);
    plot(x,wu,'--','linewidth',2)

    hold off

elseif(dims==2)

    figure; imagesc(reshape(U(end,ui),m,m)); % colorbar;
    
    clim([0 32]); hold on
    
    axis off
end