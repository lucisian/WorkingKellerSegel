clear
%set the number of dimensions

%parameters in the reaction kinematics
a=81;
D=1;
p2=1;
p4=1;
epsilon=0.3;

%load the parameters from the parameters.m script
run('parameters2.m')

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
    chi = @(U)U./(1+U.^(2));
    F = @(t,U)[f(U(ui),U(vi)) + Lap*U(ui) -...
    p*(chi(U(ui)).*(Lap*U(vi)) + (Adv*(chi(U(ui)))).*(Adv*U(vi)));
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

    A03=-6*sqrt((sqrt(a)*p2/((sqrt(a)+sqrt(D))*(9*a-67*sqrt(a*D)-16*D))));

    w0=1;
    w11=1;
    w22=(-(a+4*sqrt(a*D)))/(18*a);
    w20=-1/2;
    w=w0+epsilon*A03*w11*cos(k*x) +epsilon*epsilon*w22*A03*A03*cos(2*k*x) +epsilon*epsilon*w20*A03*A03;

    plot(x,w,'--','linewidth',2)

    hold on
    
    A0=-0.453859;
    wu=1-A0^(2)*epsilon^(2)/2 -267925*A0^(4)*epsilon^(4)/662661 +899*A0^(2)*epsilon^(5)/20200 ...
    +(A0*epsilon -55*A0^(3)*epsilon^(3)/32724 +A0*epsilon^(4)/2020)*cos(k*x) ...
    +(-13*A0^(2)*epsilon^(2)/162 +8956342085*A0^(4)*epsilon^(4)/6870469248 -79109*A0^(2)*epsilon^(5)/132532200)*cos(2*k*x) ...
    +(-1189*A0^(3)*epsilon^(3)/10368)*cos(3*k*x) ...
    +(3385445*A0^(4)*epsilon^(4)/68024448)*cos(4*k*x);
    plot(x,wu,'--','linewidth',2)

    hold off

elseif(dims==2)

    figure; imagesc(reshape(U(end,ui),m,m)); % colorbar;
    
    clim([0 32]); hold on
    
    axis off
end