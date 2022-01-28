% Extra material for the paper
% General construction and classes of explicit L1-optimal couplings
% by G. Puccetti ane L. Rüschendorf

% This code produces Figure 1.1 
% showing optimal transportation lines with respect to the squared Euclidean distance (p = 2, left) 
% and the Euclidean distance (p = 1, right), between the same Gaussian distributions 
% (corresponding to case B in Table 6.1). The solution for p = 1 is given for the first time in the paper.

% The parameters of the following code are calibrated to obtain Figure 1.1
% Code might not properly work and might require adjustements if applied to a different example

% auxiliary files needed: subtightplot.m / condition2.m/ condition3.m

clear;
%suptightplot is called (the .m file must be included in the same directory)
%for nicer plots
opt = {0.03, 0.03, 0.03} ;
subplot = @(m,n,p) subtightplot(m,n,p,opt{:}); 
figure();
%resol=number of transportation rays and resolution of the integration grid
resol=25;
%line widht
lin=1.7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LEFT PLOT - L2 case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,1);
figure(1);
%%% NUMBER OF TRANSPORTATION RAYS (the more the rays the more accurate
%%% the final estimate, the more time consuming the code)
% upper bound for the domain of a Gaussian distribution
% everything will be restricted to [-M,M]^2
M=4;
%%% INPUTS: COVARIANCE MATRICES of the two Guassian measures %%%
% mean vectors are assumed to be null
mux = [0 0];
muy = [0 0];
global xdiag SX SY sx sy d1 d2 int1 int2 delta1 delta2;
%upperbound for unbounded numerical integrals 
MMM=6;
%Covariance matrix of X,Y. They must be invertible.
SX = [1 -0.4; -0.4 1];
SY = [1 0.8; 0.8 1];
%Plot the two diagonals of the plane
NN=400;
u=linspace(-M,M,NN);
plot(u,u,'k','LineWidth',2)
hold on
plot(u,-u,'k','LineWidth',2)
xticks([]);
yticks([]);
axis square
hold on;
% trasportation lines in the L2 case
% grid
MMMM=30;
xgrid=linspace(-M,M,MMMM);
%L2-optimal transportation (see for instance Corollary 3.2.13 in Rachev and
%Rüschendorf (1998))
T=inv(sqrtm(SX))*sqrtm(sqrtm(SX)*SY*sqrtm(SX))*inv(sqrtm(SX));
%Computation of destination point for arrows
N=1;
for i=1:MMMM
    for j=1:MMMM
    xstart(N)=xgrid(i);
    ystart(N)=xgrid(j);    
    xdest(N)=T(1,1)*xstart(N)+T(1,2)*ystart(N);
    ydest(N)=T(2,1)*xstart(N)+T(2,2)*ystart(N);
    N=N+1;     
    end
end
%plot arrows
for i=1:MMMM^2
plot([xstart(i), xdest(i)], [ystart(i), ydest(i)],'k','LineWidth',lin)
axis([-M M -M M])
hold on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RIGHT PLOT - L1 case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,2);
%Plot the two diagonals of the plane
NN=400;
u=linspace(-M,M,NN);
plot(u,u,'k','LineWidth',2)
hold on
plot(u,-u,'k','LineWidth',2)
xticks([]);
yticks([]);
axis square
hold on;
%%% NUMBER OF TRANSPORTATION RAYS (the more the rays the more accurate
%%% the final estimate the more time consuming the code)

%determinants of covariance matrices
dx=det(SX);
dy=det(SY);
%compute the rotation which is necessary to 
%make A and B symmetric wrt to the sectors as in the paper
%constant as in the paper
V=(SX(2,2)-SX(1,1))/dx-(SY(2,2)-SY(1,1))/dy;
W=SY(1,2)/dy-SX(1,2)/dx;
%varrho=rotation angle (if angle is 0 no rotation is needed)
varrho=atan(V/2/W)/2;
%definition of rotation matrix
D=[cos(varrho) -sin(varrho); sin(varrho) cos(varrho)];
%new covariance matrices after rotation
SX=D*SX*transpose(D);
SY=D*SY*transpose(D);
%correlation parameters to be used in the density formulas
sx=SX(1,2)/sqrt(SX(1,1)*SX(2,2));
sy=SY(1,2)/sqrt(SY(1,1)*SY(2,2));

%%% %%% %%% %%% TRANSPORTATION RAYS AND  COST %%% %%% %%% %%% %%% %%% %%% 

% Transportation rays are identified by dividing sector S1 
% into many slices so that each slice has the same probability according to
% the two measures
 
% identification of one point on the main diagonal (x,x)
% between epsilon and M
epsilon=0.008; 
%rays starting points on the main diagonal
xt=linspace(epsilon,M-epsilon,resol);
deltaxt=xt(2)-xt(1);
%auxiliary vectors/parameters
dir1=linspace(0,0,resol);
dir2=linspace(0,0,resol);
slope=linspace(0,0,resol);
upp=linspace(0,0,resol);
ddd=linspace(0,0,resol);
cost=zeros((resol-1),(resol-1));

%bounds on the optimal direction angle
eps=0.00000001;
lowerbound=5/4*pi+eps;
upperbound=7/4*pi-eps;

for t=1:resol
    t
%starting point on the main diagonal
xdiag=xt(t);
lowerbound=5/4*pi+eps;
upperbound=7/4*pi-eps;
%direction of optimal transport computed by taking slices
rho_0=fzero(@condition2,[lowerbound,upperbound]);
dir1(t)=cos(rho_0);
dir2(t)=sin(rho_0);
upp(t)=-2*xdiag/(dir1(t)+dir2(t));
end


%%%% OPTIONAL: PLOT THE CANDIDATE TRANSPORTATION RAY AND SHOW
%%%% THAT NECESSARY CONDITION IS SATISFIED IN THE LIMIT OF RESOL->INF

s0b=linspace(0,0,(resol-1));
%calculation of boundary points
for t=1:(resol-1)
%point on the upper diagonal
xdiag=xt(t);
%recall directions
d1=dir1(t);
d2=dir2(t);
%calculation of boundary point
s0b(t)=fzero(@condition3,[0,MMM]);
%plot(xdiag+s0b(t)*d1,xdiag+s0b(t)*d2,'ko');
end
s0max=max(s0b);
figure(1);
for t=1:(resol-1)
    %print t to know code advance
    t
%point on the upper diagonal
xdiag=xt(t);
%recall directions
d1=dir1(t);
d2=dir2(t);
%discrete derivative of the two directions in t
%(mathematical details are given in the notes or in the paper)
delta1=(dir1(t+1)-dir1(t))/deltaxt;
delta2=(dir2(t+1)-dir2(t))/deltaxt;
upp=-2*xdiag/(d1+d2);
%functions for conditional densities
f1 = @(u) (abs(d2-d1+u*(delta1*d2-delta2*d1)).*(1/(2*pi*sqrt(SX(1,1)*SX(2,2)*(1-sx^2)))).*(exp((-1/(2*(1-sx^2))).*((xdiag+u.*d1).^2/SX(1,1)+(xdiag+u.*d2).^2/SX(2,2)-2.*sx.*((xdiag+u.*d1)/sqrt(SX(1,1))).*((xdiag+u.*d2)/sqrt(SX(2,2)))))));
f2 = @(u) (abs(d2-d1+u*(delta1*d2-delta2*d1)).*(1/(2*pi*sqrt(SY(1,1)*SY(2,2)*(1-sy^2)))).*(exp((-1/(2*(1-sy^2))).*((xdiag+u.*d1).^2/SY(1,1)+(xdiag+u.*d2).^2/SY(2,2)-2.*sy.*((xdiag+u.*d1)/sqrt(SY(1,1))).*((xdiag+u.*d2)/sqrt(SY(2,2)))))));
%densities of the two measures
den1 = @(u) (1/(2*pi*sqrt(SX(1,1)*SX(2,2)*(1-sx^2)))).*(exp((-1/(2*(1-sx^2))).*((xdiag+u.*d1).^2/SX(1,1)+(xdiag+u.*d2).^2/SX(2,2)-2.*sx.*((xdiag+u.*d1)/sqrt(SX(1,1))).*((xdiag+u.*d2)/sqrt(SX(2,2))))));
den2 = @(u) (1/(2*pi*sqrt(SY(1,1)*SY(2,2)*(1-sy^2)))).*(exp((-1/(2*(1-sy^2))).*((xdiag+u.*d1).^2/SY(1,1)+(xdiag+u.*d2).^2/SY(2,2)-2.*sy.*((xdiag+u.*d1)/sqrt(SY(1,1))).*((xdiag+u.*d2)/sqrt(SY(2,2))))));
%densities of the two reduced measures
den1r = @(u) max((den1(u)-den2(u)),0);
den2r = @(u) max((den2(u)-den1(u)),0);
%conditional densities for the reduced measures
g1 = @(u) max((f1(u)-f2(u)),0);
g2 = @(u) max((f2(u)-f1(u)),0);
%corresponding integrals and their difference 
int1 = integral(g1,0,upp);
int2 = integral(g2,0,upp);
%check of necessary condition; this has to go to zero as resol->Inf
ddd(t)=int1-int2;
%plot slices on the above produced plot
tseq=linspace(0,upp,resol);
xp=linspace(0,0,resol);
yp=linspace(0,0,resol);
for i=1:resol
xp(i)=xdiag+tseq(i)*d1;
yp(i)=xdiag+tseq(i)*d2;
end
hold on;
plot(xp,yp,'k','LineWidth',lin)
axis square;
hold on;
%plot slices in the other sectors
plot(yp,xp,'k','LineWidth',lin)
hold on;
plot(-yp,-xp,'k','LineWidth',lin)
hold on;
plot(-xp,-yp,'k','LineWidth',lin)
hold off;
end