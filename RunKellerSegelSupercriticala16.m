clear
%set the number of dimensions

%parameters in the reaction kinematics
a=16;
D=1;
p2=1;
p4=1;
epsilon=0.1;

%load the parameters from the parameters.m script
run('parameters.m')

% (Sparse) Laplacian matrix
e=ones(m,1);
Lap=spdiags([e,-2*e,e],[1,0,-1],m,m);
Adv=spdiags([e,-e],[1,-1],m,m);

%periodic boundary conditions
%Lap(1,end)=1;
%Lap(end,1)=1;
%Adv(1,end)=-1;
%Adv(end,1)=1;

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

%TODO: WRITE THE BELOW FOR A GENERAL SENSITIVITY FUNCTION AND BE SURE WHAT
%SENSITIVITY YOU ARE USING!!!!!!

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

%Initial condition - this is a small normally distributed perturbation
% of the homogeneous steady state of our kinetics
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

    A03=-6*sqrt((sqrt(a*D)*p2)/((sqrt(a)+sqrt(D))*(-2*a*sqrt(a)+23*a*sqrt(D)+5*sqrt(a)*D+10*sqrt(D)*D)));

    w0=1;
    w11=1;
    w22=(2*sqrt(a)+sqrt(D)*(sqrt(a)+4*sqrt(D)))/(18*sqrt(a*D));
    w20=-1/2;
    w=w0+epsilon*A03*w11*cos(k*x) +epsilon*epsilon*w22*A03*A03*cos(2*k*x) +epsilon*epsilon*w20*A03*A03;

    plot(x,w,'--','linewidth',2)

    hold on
    
    A0=-0.328148;
    wu=1-A0^(2)*epsilon^(2)/2 +99*A0^(2)*epsilon^(4)/650 -453*A0^(4)*epsilon^(4)/208 ...
        +(A0*epsilon +A0*epsilon^(3)/130 -15*A0^(3)*epsilon^(3)/208)*cos(k*x)...
        +(A0^(2)*epsilon^(2) +343*A0^(2)*epsilon^(4)/11700 -979*A0^(4)*epsilon^(4)/39936)*cos(2*k*x)...
        +(1651*A0^(3)*epsilon^(3)/2048)*cos(3*k*x)...
        +(3899*A0^(4)*epsilon^(4)/6144)*cos(4*k*x);
    plot(x,wu,'--','linewidth',2)

    hold off

elseif(dims==2)

    figure; imagesc(reshape(U(end,ui),m,m)); % colorbar;
    
    clim([0 32]); hold on
    
    axis off
end