clear

%parameters in the reaction kinematics
a=4;
D=1;
p2=1;
p4=0;
epsilon=0.2;

%load the parameters from the parameters.m script
run('parameters4.m')

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

if(dims==1)
    plot(x,U(end,ui),'linewidth',3); hold on

    set(gca,'fontsize',24);

elseif(dims==2)

    figure; imagesc(reshape(U(end,ui),m,m)); % colorbar;
    
    clim([0 32]); hold on
    
    axis off
end

[X, t] = meshgrid(x, T);

% Create the plot
figure;
imagesc(T, x, U(1:end,1:end/2)');
colorbar;
xlabel('Time');
ylabel('Space');

% Set up for frequency tracking
fs = m; % sampling frequency in space
Nl = m; % number of spatial points
freq = 0:fs/Nl:fs/2; % frequency axis (positive freqs only)
dominant_freq = zeros(length(T),1); % to store dominant frequency at each time

% Loop through each time step
for i = 1:length(T)
    u = U(i, ui); % extract u(x) at current time
    u = u - mean(u); % zero-mean for FFT
    
    xdft = fft(u);
    xdft = xdft(1:Nl/2+1); % one-sided FFT
    psdx = (1/(fs*Nl)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    
    % Find dominant frequency (highest power)
    [~, idx] = max(psdx);
    dominant_freq(i) = freq(idx);
end

% Plot the dominant frequency over time
figure;
plot(T, 2*dominant_freq, 'LineWidth', 2);
xlabel('Time');
ylabel('Dominant Frequency');
title('Dominant Spatial Frequency Over Time');
grid on;