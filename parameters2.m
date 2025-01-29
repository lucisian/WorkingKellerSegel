%parameters2.m
%Define the parameters that we do not frequently very

%Set the number of dimensions (1D or 2D)
dims=1;

%Number of gridpoints per dimension
m=4000;

%Total number of gridpoints; varies by dimension as:
%2D make N=m^2; 1D make N=m;
if(dims==1)
N=m;
elseif(dims==2)
N=m^2;
end

%parameters in the reaction kinematics
p0=2*(a+D+2*sqrt(a*D));
p=p0+epsilon*epsilon*p2+epsilon*epsilon*epsilon*epsilon*p4;

%Domain length
L=3*2*pi*sqrt(sqrt(D/a));

%Spatial step size
dx=L/(m-1);

%Spatial domain (needed for plotting only)
x=linspace(0,L,m);

%Time interval to solve the equations on
T=linspace(0,10000,100);

k=sqrt(sqrt(a/D));