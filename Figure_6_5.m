% Extra material for the paper
% General construction and classes of explicit L1-optimal couplings
% by G. Puccetti ane L. Rüschendorf

% This code produces Figure 6.5  which illustrates transportation lines
% in six different cases as in Table 6.1

% this code is actually a copy and merge of the single codes that one
% can retrieve for the spcific cases.

% auxiliary files needed: subtightplot.m

clear;
format long
%suptightplot is called (must be included in the same directory)
%for nicer plots
opt = {0.03, 0.03, 0.03} ;
subplot = @(m,n,p) subtightplot(m,n,p,opt{:}); 
figure();

%resol=number of transportation rays and resolution of the integration grid
resol=15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOT 1 OF 6 - symmetric case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,1);


%%% NUMBER OF TRANSPORTATION RAYS (the more the rays the more accurate
%%% the final estimate the more time consuming the code)
% upper bound for the domain of a Gaussian distribution
% everything will be restricted to [-M,M]^2
M=4;
%%% INPUTS: COVARIANCE MATRICES of the two Gaussian measures %%%
%mean vectors are assumed to be null
mux = [0 0];
muy = [0 0];
global xdiag SX SY sx sy d1 d2 int1 int2 s1 s0 delta1 delta2;
%upperbound for integrals 
MMM=6;
%covariance matrix of X,Y. They must be invertible.
SX = [1 -0.4; -0.4 1];
SY = [1 0.4; 0.4 1];
%determinants of covariance matrices
dx=det(SX);
dy=det(SY);
%check that one is in the unbounded domain case
if prod(eig(inv(SY)-inv(SX)))>=0
error('Covariance matrices correspond to the bounded case. Please use the other code.');
      return;
end
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

%%% %%% %%% PLOT OF THE DOMAIN OF THE REDUCED DENSITIES %%% %%% %%% %%% 
NN=400;
u=linspace(-M,M,NN);
v=linspace(-M,M,NN);
fun=zeros([NN NN]);
%ddif=difference between the densities of the two original distributions
ddif= @(x,y) max(0,(mvnpdf([x,y],[],SX)-mvnpdf([x,y],[],SY)));
%assign different colors where open density is bigger/smaller than the other
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
%plot the two diagonals
hold on
plot(u,u,'k','LineWidth',2)
hold on
plot(u,-u,'k','LineWidth',2)
xticks([]);
yticks([]);
axis square
hold on;
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%


%%% %%% %%% %%% TRANSPORTATION RAYS AND  COST %%% %%% %%% %%% %%% %%% %%% 

% transportation rays are identified by dividing sector S1 
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

for t=1:resol    t
%starting point on the main diagonal
xdiag=xt(t);
lowerbound=5/4*pi+eps;
upperbound=7/4*pi-eps;
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
plot(xp,yp,'k','LineWidth',2)
axis square;
hold on;
%plot slices in the other sectors
plot(yp,xp,'k','LineWidth',2)
hold on;
plot(-yp,-xp,'k','LineWidth',2)
hold on;
plot(-xp,-yp,'k','LineWidth',2)
hold on;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOT 2 OF 6 - unbounded case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,2);

%%% NUMBER OF TRANSPORTATION RAYS (the more the rays the more accurate
%%% the final estimate the more time consuming the code)
% upper bound for the domain of a Gaussian distribution
% everything will be restricted to [-M,M]^2
M=4;
%%% INPUTS: COVARIANCE MATRICES of the two Guassian measures %%%
%mean vectors are assumed to be null
mux = [0 0];
muy = [0 0];
global xdiag SX SY sx sy d1 d2 int1 int2 s1 s0 delta1 delta2;
%upperbound for integrals 
MMM=6;
%covariance matrix of X,Y. They must be invertible.
SX = [1 -0.4; -0.4 1];
SY = [1 0.8; 0.8 1];
%determinants of covariance matrices
dx=det(SX);
dy=det(SY);
%check that one is in the unbounded domain case
if prod(eig(inv(SY)-inv(SX)))>=0
error('Covariance matrices correspond to the bounded case. Please use the other code.');
      return;
end
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

%%% %%% %%% PLOT OF THE DOMAIN OF THE REDUCED DENSITIES %%% %%% %%% %%% 
NN=400;
u=linspace(-M,M,NN);
v=linspace(-M,M,NN);
fun=zeros([NN NN]);
%ddif=difference between the densities of the two original distributions
ddif= @(x,y) max(0,(mvnpdf([x,y],[],SX)-mvnpdf([x,y],[],SY)));
%assign different colors where open density is bigger/smaller than the other
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
% also plot the two diagonals
hold on
plot(u,u,'k','LineWidth',2)
hold on
plot(u,-u,'k','LineWidth',2)
xticks([]);
yticks([]);
axis square
hold on;
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%


%%% %%% %%% %%% TRANSPORTATION RAYS AND  COST %%% %%% %%% %%% %%% %%% %%% 

% transportation rays are identified by dividing sector S1 
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
plot(xp,yp,'k','LineWidth',2)
axis square;
hold on;
%plot slices in the other sectors
plot(yp,xp,'k','LineWidth',2)
hold on;
plot(-yp,-xp,'k','LineWidth',2)
hold on;
plot(-xp,-yp,'k','LineWidth',2)
hold on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOT 3 OF 6 - unbounded case 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,3);


%%% NUMBER OF TRANSPORTATION RAYS (the more the rays the more accurate
%%% the final estimate the more time consuming the code)
% upper bound for the domain of a Gaussian distribution
% everything will be restricted to [-M,M]^2
M=4;
%%% INPUTS: COVARIANCE MATRICES of the two Guassian measures %%%
%mean vectors are assumed to be null
mux = [0 0];
muy = [0 0];
global xdiag SX SY sx sy d1 d2 int1 int2 s1 s0 delta1 delta2;
%upperbound for integrals 
MMM=6;
%covariance matrix of X,Y. They must be invertible.
%1.8 is approximated by 1.799 otherwise use Gaussian_OC_caseC
SX = [1.799 -0.4; -0.4 1.799];
SY = [1 0.4; 0.4 1];
%determinants of covariance matrices
dx=det(SX);
dy=det(SY);
%check that one is in the unbounded domain case
if prod(eig(inv(SY)-inv(SX)))>=0
error('Covariance matrices correspond to the bounded case. Please use the other code.');
      return;
end
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

%%% %%% %%% PLOT OF THE DOMAIN OF THE REDUCED DENSITIES %%% %%% %%% %%% 
NN=400;
u=linspace(-M,M,NN);
v=linspace(-M,M,NN);
fun=zeros([NN NN]);
%ddif=difference between the densities of the two original distributions
ddif= @(x,y) max(0,(mvnpdf([x,y],[],SX)-mvnpdf([x,y],[],SY)));
%assign different colors where open density is bigger/smaller than the other
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
% also plot the two diagonals
hold on
plot(u,u,'k','LineWidth',2)
hold on
plot(u,-u,'k','LineWidth',2)
xticks([]);
yticks([]);
axis square
hold on;
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%


%%% %%% %%% %%% TRANSPORTATION RAYS AND  COST %%% %%% %%% %%% %%% %%% %%% 

% transportation rays are identified by dividing sector S1 
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
plot(xp,yp,'k','LineWidth',2)
axis square;
hold on;
%plot slices in the other sectors
plot(yp,xp,'k','LineWidth',2)
hold on;
plot(-yp,-xp,'k','LineWidth',2)
hold on;
plot(-xp,-yp,'k','LineWidth',2)
hold on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOT 4 OF 6 - bounded case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,4);
% upper bound for the domain of a Gaussian distribution
% everything will be restricted to [-M,M]^2
M=5;
%%% INPUTS: COVARIANCE MATRICES of the two Guassian measures %%%
%mean vectors are assumed to be null
mux = [0 0];
muy = [0 0];
global xdiag SX SY sx sy d1 d2 int1 int2 s1 s0 delta1 delta2 MMM xell yell;
%%%Upper bound on unbounded integrals
MMM=10;
%covariance matrix of X,Y. They must be invertible.
SX = [2 -0.4; -0.4 2];
SY = [1 0.4; 0.4 1];
%determinants of covariance matrices
dx=det(SX);
dy=det(SY);
%check that one is in the unbounded domain case
if prod(eig(inv(SY)-inv(SX)))<=0
error('Covariance matrices correspond to the UNbounded case. Please use the other code.');
      return;
end
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
% % % % % % % % % % % % %correlation parameters to be used in the density formulas
sx=SX(1,2)/sqrt(SX(1,1)*SX(2,2));
sy=SY(1,2)/sqrt(SY(1,1)*SY(2,2));

%%% %%% %%% PLOT OF THE DOMAIN OF THE REDUCED DENSITIES %%% %%% %%% %%% 
NN=400;
u=linspace(-M,M,NN);
v=linspace(-M,M,NN);
fun=zeros([NN NN]);
%ddif=difference between the densities of the two original distributions
ddif= @(x,y) max(0,(mvnpdf([x,y],[],SX)-mvnpdf([x,y],[],SY)));
%assign different colors where open density is bigger/smaller than the other
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


%%% %%% %%% %%% TRANSPORTATION RAYS AND  COST %%% %%% %%% %%% %%% %%% %%% 

% transportation rays are identified by dividing sector S1 
% into many slices so that each slice has the same probability according to
% the two measures

% wrt to the unbounded case, we discretize the ellipse represinting the boundary between reduced
% densities domains. We use a number of rays equal to resol

%find a discretization over the ellipse
angle=linspace(7/4*pi,9/4*pi, resol+2);
angle(1)=[];
angle(resol+1)=[];
elpoint=linspace(0,0,resol);

%auxiliary viarables
%optimal directions
dir1=linspace(0,0,resol);
dir2=linspace(0,0,resol);
xt=linspace(0,0,resol);
cost=zeros(resol-1,resol-1);


%for each angle find boundary point on the ellipse
%eps=tolerance parameter for uniroot at the boundaries
eps=0.0001;
%here loop on i 
for i=1:resol
xdiag=0;
d1=cos(angle(i));
d2=sin(angle(i));
%find boundary point on the ellipse
elpoint(i)=fzero(@condition3,[0,MMM]);
%coordinate of point on the ellipse
xell=xdiag+elpoint(i)*d1;
yell=xdiag+elpoint(i)*d2;
%plot(xell,yell,'ro');
%find optimal directions
rho0=fzero(@condition2b,[7/4*pi+eps,9/4*pi-eps]);
dir1(i)=cos(rho0);
dir2(i)=sin(rho0);
%corresponding points on the diagonal to form first dimension of
%integration grid
xt(i)=xell+dir1(i)*(xell-yell)/(dir2(i)-dir1(i));
%plot slices on the above produced plot
tseq=linspace(0,MMM,resol);
xp=linspace(0,0,resol);
yp=linspace(0,0,resol);
for j=1:resol
xp(j)=xt(i)+tseq(j)*dir1(i);
yp(j)=xt(i)+tseq(j)*dir2(i);
end
figure(1)
hold on;
plot(xp,yp,'k','LineWidth',2)
axis square;
hold on;
%plot slices in the other sectors
plot(yp,xp,'k','LineWidth',2)
hold on;
plot(-yp,-xp,'k','LineWidth',2)
hold on;
plot(-xp,-yp,'k','LineWidth',2)
hold on;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOT 5 OF 6 - bounded case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,5);

% upper bound for the domain of a Gaussian distribution
% everything will be restricted to [-M,M]^2
M=5;
%%% INPUTS: COVARIANCE MATRICES of the two Guassian measures %%%
%mean vectors are assumed to be null
mux = [0 0];
muy = [0 0];
global xdiag SX SY sx sy d1 d2 int1 int2 s1 s0 delta1 delta2 MMM xell yell;
%%%Upper bound on unbounded integrals
MMM=10;
%covariance matrix of X,Y. They must be invertible.
SX = [1.3 0.4; 0.4 1.3];
SY = [1 0.4; 0.4 1];
%determinants of covariance matrices
dx=det(SX);
dy=det(SY);
%check that one is in the unbounded domain case
if prod(eig(inv(SY)-inv(SX)))<=0
error('Covariance matrices correspond to the UNbounded case. Please use the other code.');
      return;
end
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
% % % % % % % % % % % % %correlation parameters to be used in the density formulas
sx=SX(1,2)/sqrt(SX(1,1)*SX(2,2));
sy=SY(1,2)/sqrt(SY(1,1)*SY(2,2));

%%% %%% %%% PLOT OF THE DOMAIN OF THE REDUCED DENSITIES %%% %%% %%% %%% 
NN=400;
u=linspace(-M,M,NN);
v=linspace(-M,M,NN);
fun=zeros([NN NN]);
%ddif=difference between the densities of the two original distributions
ddif= @(x,y) max(0,(mvnpdf([x,y],[],SX)-mvnpdf([x,y],[],SY)));
%fun3 = @(x,y) max(0,(((1/(2*pi*sqrt(1-sx^2))).*(exp((-1/(2*(1-sx^2))).*(x.^2+y.^2-2.*sx.*x.*y))))-((1/(2*pi*sqrt(1-sy^2))).*(exp((-1/(2*(1-sy^2))).*(x.^2+y.^2-2.*sy.*x.*y))))));
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
% also plot the two diagonals
hold on
plot(u,u,'k','LineWidth',2)
hold on
plot(u,-u,'k','LineWidth',2)
xticks([]);
yticks([]);
axis square
hold on;
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%


%%% %%% %%% %%% TRANSPORTATION RAYS AND  COST %%% %%% %%% %%% %%% %%% %%% 

% transportation rays are identified by dividing sector S1 
% into many slices so that each slice has the same probability according to
% the two measures

% wrt to the unbounded case, we discretize the ellipse represinting the boundary between reduced
% densities domains. We use a number of rays equal to resol

%find a discretization over the ellipse
angle=linspace(7/4*pi,9/4*pi, resol+2);
angle(1)=[];
angle(resol+1)=[];
elpoint=linspace(0,0,resol);

%auxiliary viarables
%optimal directions
dir1=linspace(0,0,resol);
dir2=linspace(0,0,resol);
xt=linspace(0,0,resol);
cost=zeros(resol-1,resol-1);


%for each angle find boundary point on the ellipse
%eps=tolerance parameter for uniroot at the boundaries
eps=0.0001;
%here loop on i 
for i=1:resol
xdiag=0;
d1=cos(angle(i));
d2=sin(angle(i));
%find boundary point on the ellipse
elpoint(i)=fzero(@condition3,[0,MMM]);
%coordinate of point on the ellipse
xell=xdiag+elpoint(i)*d1;
yell=xdiag+elpoint(i)*d2;
%plot(xell,yell,'ro');
%find optimal directions
rho0=fzero(@condition2b,[7/4*pi+eps,9/4*pi-eps]);
dir1(i)=cos(rho0);
dir2(i)=sin(rho0);
%corresponding points on the diagonal to form first dimension of
%integration grid
xt(i)=xell+dir1(i)*(xell-yell)/(dir2(i)-dir1(i));
%plot slices on the above produced plot
tseq=linspace(0,MMM,resol);
xp=linspace(0,0,resol);
yp=linspace(0,0,resol);
for j=1:resol
xp(j)=xt(i)+tseq(j)*dir1(i);
yp(j)=xt(i)+tseq(j)*dir2(i);
end
figure(1)
hold on;
plot(xp,yp,'k','LineWidth',2)
axis square;
hold on;
%plot slices in the other sectors
plot(yp,xp,'k','LineWidth',2)
hold on;
plot(-yp,-xp,'k','LineWidth',2)
hold on;
plot(-xp,-yp,'k','LineWidth',2)
hold on;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOT 6 OF 6 - radial case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,6);

% upper bound for the domain of a Gaussian distribution
% everything will be restricted to [-M,M]^2
M=5;
%%% INPUTS: COVARIANCE MATRICES of the two Guassian measures %%%
%mean vectors are assumed to be null
mux = [0 0];
muy = [0 0];
global xdiag SX SY sx sy d1 d2 int1 int2 s1 s0 delta1 delta2 MMM xell yell;
%%%Upper bound on unbounded integrals
MMM=10;
%covariance matrix of X,Y. They must be invertible.
SX = [2 0; 0 2];
SY = [1 0; 0 1];
% % % % % % % % % % % % %correlation parameters to be used in the density formulas
sx=SX(1,2)/sqrt(SX(1,1)*SX(2,2));
sy=SY(1,2)/sqrt(SY(1,1)*SY(2,2));

%%% %%% %%% PLOT OF THE DOMAIN OF THE REDUCED DENSITIES %%% %%% %%% %%% 
NN=400;
u=linspace(-M,M,NN);
v=linspace(-M,M,NN);
fun=zeros([NN NN]);
%ddif=difference between the densities of the two original distributions
ddif= @(x,y) max(0,(mvnpdf([x,y],[],SX)-mvnpdf([x,y],[],SY)));
%assign different colors where open density is bigger/smaller than the other
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
% also plot the two diagonals
hold on
plot(u,u,'k','LineWidth',2)
hold on
plot(u,-u,'k','LineWidth',2)
xticks([]);
yticks([]);
axis square
hold on;
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%


%%% %%% %%% %%% TRANSPORTATION RAYS AND  COST %%% %%% %%% %%% %%% %%% %%% 

% transportation rays are identified by dividing sector S1 
% into many slices so that each slice has the same probability according to
% the two measures

% wrt to the unbounded case, we discretize the ellipse represinting the boundary between reduced
% densities domains. We use a number of rays equal to resol

%find a discretization over the ellipse
angle=linspace(7/4*pi,9/4*pi, resol+2);
angle(1)=[];
angle(resol+1)=[];
elpoint=linspace(0,0,resol);

%auxiliary viarables
%optimal directions
dir1=linspace(0,0,resol);
dir2=linspace(0,0,resol);
xt=linspace(0,0,resol);
cost=zeros(resol-1,resol-1);


%for each angle find boundary point on the ellipse
%eps=tolerance parameter for uniroot at the boundaries
eps=0.0001;
%here loop on i 
for i=1:resol
xdiag=0;
d1=cos(angle(i));
d2=sin(angle(i));
%find boundary point on the ellipse
elpoint(i)=fzero(@condition3,[0,MMM]);
%coordinate of point on the ellipse
xell=xdiag+elpoint(i)*d1;
yell=xdiag+elpoint(i)*d2;
%plot(xell,yell,'ro');
%find optimal directions
rho0=fzero(@condition2b,[7/4*pi+eps,9/4*pi-eps]);
dir1(i)=cos(rho0);
dir2(i)=sin(rho0);
%corresponding points on the diagonal to form first dimension of
%integration grid
xt(i)=xell+dir1(i)*(xell-yell)/(dir2(i)-dir1(i));
%plot slices on the above produced plot
tseq=linspace(0,MMM,resol);
xp=linspace(0,0,resol);
yp=linspace(0,0,resol);
for j=1:resol
xp(j)=xt(i)+tseq(j)*dir1(i);
yp(j)=xt(i)+tseq(j)*dir2(i);
end
figure(1)
hold on;
plot(xp,yp,'k','LineWidth',2)
axis square;
hold on;
%plot slices in the other sectors
plot(yp,xp,'k','LineWidth',2)
hold on;
plot(-yp,-xp,'k','LineWidth',2)
hold on;
plot(-xp,-yp,'k','LineWidth',2)
hold on;
end

hold off;
subplot(3,2,1);
title(['A (symmetric case)'],'FontSize', 14);
subplot(3,2,2);
title(['B (unbounded case)'],'FontSize', 14)
subplot(3,2,3);
title(['C (unbounded case)'],'FontSize', 14);
subplot(3,2,4);
title(['D (bounded case)'],'FontSize', 14)
subplot(3,2,5);
title(['E (bounded case)'],'FontSize', 14)
subplot(3,2,6);
title(['F (radial case)'],'FontSize', 14)


hold off;