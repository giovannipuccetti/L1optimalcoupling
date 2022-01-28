% Extra material for the paper
% General construction and classes of explicit L1-optimal couplings
% by G. Puccetti ane L. Rüschendorf

% This code approximates the trasportation cost given by the coupling described
% in the paper in the case B in Table 6.1,

% The parameters of the following code are calibrated to the specific case
% Code might not properly work and might require adjustements if applied to a different example

% auxiliary files needed: condition2.m /condition3.m /condition4.m

%this version seems to be obsolete and replaced by new_unbounded_case

%results
%0.748175684240197 with M=4 and resol=500 points

clear;

%%% NUMBER OF TRANSPORTATION RAYS (the more the rays the more accurate
%%% the final estimate the more time consuming the code)

%resol=number of transportation rays and resolution of the integration grid
resol=500;
% Gaussian distribution will be restricted to [0,M]^2
M=4;
%%% INPUTS: COVARIANCE MATRICES of the two Gaussian measures %%%
%   (mean vectors are assumed to be null)
mux = [0 0];
muy = [0 0];
global xdiag SX SY sx sy d1 d2 int1 int2 s1 s0 delta1 delta2;
%upperbound for numerical integrals 
MMM=10;
%covariance matrix of X,Y. They must be invertible.
SX = [1 -0.4; -0.4 1];
SY = [1 0.8; 0.8 1];
%determinants of covariance matrices
dx=det(SX);
dy=det(SY);
%check that the two covariance matrices imply reduced measures with
%unbounded domains
if prod(eig(inv(SY)-inv(SX)))>=0
error('Covariance matrices correspond to the bounded case. Please use the other code Gaussain_ bounded_case.');
      return;
end
%Compute the rotation which is necessary to 
%make the domains A and B of the reduced densities symmetric wrt 
%to the sectors 
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

%%% %%% %%% PLOT OF THE DOMAIN OF THE REDUCED DENSITIES %%% %%% %%% %%% 
%grid parameter
NN=400;
u=linspace(-M,M,NN);
v=linspace(-M,M,NN);
fun=zeros([NN NN]);
%ddif=difference between the densities of the two original distributions
ddif= @(x,y) max(0,(mvnpdf([x,y],[],SX)-mvnpdf([x,y],[],SY)));
for i=1:NN
    for j=1:NN
    x=u(i);
    y=v(j);
    fun(i,j) =ddif(x,y);
    if fun(i,j)>0
        fun(i,j)=200;
    end
    if fun(i,j)==0
        fun(i,j)=500;
    end
    end
end
figure(1);
plot(30:30, 'k.-','MarkerSize',5), axis([-M M -M M]);
hold on
image([-M,M],[-M,M],fun)
colormap gray
% plot the two diagonals
hold on
plot(u,u,'k','LineWidth',2)
hold on
plot(u,-u,'k','LineWidth',2)
xticks([]);
yticks([]);
axis square
hold on;
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%

%%% %%% %%% %%% TRANSPORTATION RAYS AND TOTAL COST %%% %%% %%% %%% %%% %%% %%% 

% Optimal transportation rays are identified by dividing sector S1 
% into many slices so that each slice has the same probability according to
% the two (reduced) measures
 
% identification of one point on the main diagonal (x,x)
% between epsilon and M (choosing epsilon=0 might cause singularities)
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

%bounds on the optimal direction angle of the ray starting from (x,x)
eps=0.00000001;
lowerbound=5/4*pi+eps;
upperbound=7/4*pi-eps;
%first we find all optimal directions corresponding to each slice (and each
%xdiag)

% for t=11:11
%     xdiag=xt(t);
% xtt=linspace(lowerbound,upperbound,100);
% xzt=linspace(0,0,100);
% yzt=linspace(0,0,100);
% for i=1:100
%     yzt(i)=condition2(xtt(i));
% end
% plot(xtt,yzt,xtt,xzt);
% hold on;
% end

for t=1:resol
    t
%starting point on the main diagonal
xdiag=xt(t);
%lowerbound=5/4*pi+eps;
%upperbound=7/4*pi-eps;
%we check whether root is between lowerbound and upperbound
%if condition2(lowerbound)*condition2(upperbound)<0

%computation of direction of optimal transport
rho_0=fzero(@condition2,[lowerbound,upperbound]);
%else
    %eps=0.0000001
    %lowerbound=5/4*pi+eps;
%upperbound=7/4*pi-eps;
%rho_0=fzero(@condition2,[lowerbound,upperbound]);
%otherwise we change the bounds
% lowerbound=7/4*pi+eps;
% upperbound=9/4*pi-eps;
% rho_0=fzero(@condition2,[lowerbound,upperbound]);
%rho_0=fzero(@condition2,7/4*pi+eps);
%end
dir1(t)=cos(rho_0);
dir2(t)=sin(rho_0);
upp(t)=-2*xdiag/(dir1(t)+dir2(t));
%compute new lowerbound for next iteration
%y=xdiag-2*xdiag*dir1(t)/(dir1(t)+dir2(t));
%ang=atan((-y-xt(t+1))/(y-xt(t+1)));
%lowerbound=(ang>0)*(ang+pi)+(ang<0)*(ang+2*pi)
end
%max values of all upp (for defintiion of  the second coordinate of trapz grid)
%mupp=max(upp);

%%%% OPTIONAL: PLOT OPTIMAL TRANSPORTATION RAY AND SHOW
%%%% THAT NECESSARY CONDITION IS SATISFIED IN THE LIMIT OF RESOL->INF
s0b=linspace(0,0,(resol-1));
%calculation of s_0 boundary points for each point chosen on the diagonal
for t=1:(resol-1)
%point on the upper diagonal
xdiag=xt(t);
%recall directions
d1=dir1(t);
d2=dir2(t);
%calculation of boundary point s_0
s0b(t)=fzero(@condition3,[0,MMM]);
plot(xdiag+s0b(t)*d1,xdiag+s0b(t)*d2,'ko');
end
s0max=max(s0b);
figure(1);
for t=1:(resol-1)
%print t to know code advance
t
%point on the upper diagonal
xdiag=xt(t);
%recall optimal directions previously computed
d1=dir1(t);
d2=dir2(t);
%discrete derivative of the two directions in t
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
plot(xp,yp,'k','LineWidth',2)
axis square;
hold on;
vseq=linspace(0,s0max+1,(resol-1));
%recall s0 value
s0=s0b(t);
%computation of transportation cost for each value of the vector vseq
for tt=1:(resol-1)
%fixed s1 on such grid
s1=vseq(tt);
if s1<s0 
%if condition is satisfied s1 is moved to s2 by quantile coupling OW it is
%not moved
%xp=xdiag+vseq(tt)*d1;
%yp=xdiag+vseq(tt)*d2;
%plot(xp,yp,'ko','LineWidth',2);
s2=fzero(@condition4,[s0,upp]);
%associated cost is abs(s2-s1)
%point under consideration
%coordx(tt,t)=xdiag+s1*d1;
%coordy(tt,t)=xdiag+s1*d2;
%cost of the transport multiplied by reduced density
cost(tt,t)=abs(s2-s1)*g2(s1);
end
end
end

hold off;
figure;
%OPTIIONAL: check that necessary condition is going to zero when resol->Inf
plot(xt,ddd); 
%drop last coordinate of xt which is not used in the cost function
xt(end)=[];
%trapezoidal integration of the cost over the grid provides the final
%estimate of cost of quantile coupling multiplied by four (symmetric sectors)
format long
4*trapz(vseq,trapz(xt,cost,2))
