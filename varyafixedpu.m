clear

%parameters in the reaction kinematics
a=140;
D=1;
p2=1;
p4=0;
epsilon=0.1;

%load the parameters from the parameters.m script
run('parameters.m')

% (Sparse) Laplacian matrix
e=ones(m,1);
Lap=spdiags([e,-2*e,e],[1,0,-1],m,m);
Adv=spdiags([e,-e],[1,-1],m,m);

%Neumann
Lap(1,1)=-1;
Lap(end,end)=-1;
Adv(1,1)=-1;
Adv(end,end)=1;

%Indices corresponding to u variable and v variable.
%This lets us stack them both in a vector U and write u = U(ui) and v = U(vi).
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

%Initial condition - this is a small normally distributed perturbation
% of the homogeneous steady state of our kinetics
U0 = [1*ones(N,1);(1/a)*ones(N,1)]+1e-2*randn(2*N,1);
JacSparse=sparse([Lap,Lap;eye(N),Lap]);
odeOptions=odeset('JPattern',JacSparse);
[T,U]=ode15s(F,T,U0,odeOptions);
P=U(end,ui);
close all;

%Third order coefficient - supercritical when positive, subcritical when
%engative
C3 = (9*a-67*sqrt(a*D)-16*D)/72;

if(dims==1)
    plot(x,U(end,ui),'linewidth',2);
    axis tight
    hold on

    set(gca,'fontsize',24);

    hold on
    
    k=sqrt(sqrt(a/D));
    A0=-4.32;
    %commented out are the expressions for wu with different values of
    %parameter a

    %64wu=1 - 0.5*A0^(2)*epsilon^(2) -2.76942*A0^(4)*epsilon^(4) +0.0974104*A0^(2)*epsilon^(4)*p2...
    %   +(A0*epsilon -0.0210874*A0^(3)*epsilon^(3) +0.00135501*A0*epsilon^(3)*p2)*cos(k*x)...
    %   +(1.41667*A0^(2)*epsilon^(2) -0.693991*A0^(4)*epsilon^(4) +0.0230122*A0^(2)*epsilon^(4)*p2)*cos(2*k*x)...
    %   +1.54049*A0^(3)*epsilon^(3)*cos(3*k*x) +1.64549*A0^(4)*epsilon^(4)*cos(4*k*x);
    
    %wu=1- 0.5*A0^(2)*epsilon^(2) +0.0718139*A0^(2)*epsilon^(4) -1.91483*A0^(4)*epsilon^(4)...
    %  +(A0*epsilon +0.000479818*A0*epsilon^(3) -0.0000068095*A0^(3)*epsilon^(3))*cos(k*x)...
    %  +(1.82418*A0^(2)*epsilon^(2) +0.0227873*A0^(2)*epsilon^(4) -1.74744*A0^(4)*epsilon^(4))*cos(2*k*x)...
    %  +2.47055*A0^(3)*epsilon^(3)*cos(3*k*x) +3.27845*A0^(4)*epsilon^(4)*cos(4*k*x);

    %136wu= 1 -0.5*A0^(2)*epsilon^(2) +0.0722501*A0^(2)*epsilon^(4) -1.94968*A0^(4)*epsilon^(4)...
    %    +(A0*epsilon +0.000489556*A0*epsilon^(3) -0.000358368*A0^(3)*epsilon^(3))*cos(k*x)...
    %    +(1.81482*A0^(2)*epsilon^(2) +0.0227847*A0^(2)*epsilon^(4) -1.71984*A0^(4)*epsilon^(4))*cos(2*k*x)...
    %    +2.44684*A0^(3)*epsilon^(3)*cos(3*k*x)...
    %    +3.23271*A0^(4)*epsilon^(4)*cos(4*k*x);

    %138wu=1 -0.5*A0^(2)*epsilon^(2) +0.0713858*A0^(2)*epsilon^(4) -1.87952*A0^(4)*epsilon^(4)...
    %    +(A0*epsilon +0.000470404*A0*epsilon^(3) +0.000337895*A0^(3)*epsilon^(3))*cos(k*x)...
    %    +(1.83347*A0^(2)*epsilon^(2) +0.0227901*A0^(2)*epsilon^(4) -1.77498*A0^(4)*epsilon^(4))*cos(2*k*x)...
    %    +2.49419*A0^(3)*epsilon^(3)*cos(3*k*x) +3.32427*A0^(4)*epsilon^(4)*cos(4*k*x);

    %140
    wu=1 -0.5*A0^(2)*epsilon^(2) -1.60101*A0^(4)*epsilon^(4) +0.0684133*A0^(2)*epsilon^(4)*p2...
       +(A0*epsilon +0.00272572*A0^(3)*epsilon^(3) +0.000408742*A0*epsilon^(3)*p2)*cos(k*x)...
       +(1.90117*A0^(2)*epsilon^(2) -1.97978*A0^(4)*epsilon^(4) +0.0228139*A0^(2)*epsilon^(4)*p2)*cos(2*k*x)...
       +2.66985*A0^(3)*epsilon^(3)*cos(3*k*x) +3.6708*A0^(4)*epsilon^(4)*cos(4*k*x);

    %170wu=1 -0.5*A0^(2)*epsilon^(2) -1.30072*A0^(4)*epsilon^(4) +0.0657994*A0^(2)*epsilon^(4)*p2...
    %    +(A0*epsilon +0.00481712*A0^(3)*epsilon^(3) +0.000359624*A0*epsilon^(3)*p2)*cos(k*x)...
    %    +(1.96576*A0^(2)*epsilon^(2) -2.1814*A0^(4)*epsilon^(4) +0.0228412*A0^(2)*epsilon^(4)*p2)*cos(2*k*x)...
    %    +2.84282*A0^(3)*epsilon^(3)*cos(3*k*x) +4.02227*A0^(4)*epsilon^(4)*cos(4*k*x);


    %plot(x,wu,'--','linewidth',2)
    %axis([0 L 0.65 1.75])
    axis tight

    hold off

elseif(dims==2)

    figure; imagesc(reshape(U(end,ui),m,m)); % colorbar;
    
    clim([0 32]); hold on
    
    axis off
end