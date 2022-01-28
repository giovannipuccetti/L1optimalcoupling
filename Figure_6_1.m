% Extra material for the paper
% General construction and classes of explicit L1-optimal couplings
% by G. Puccetti ane L. Rüschendorf

% This code produces Figure 6.1 showing how the domains of two reduced Gaussian densities 
% can always be rotated (without changing their L1-distance) in order to fulfill 
% the symmetry condition as given in (6.3)

% The parameters of the following code are calibrated to obtain Figure 6.1
% Code might not properly work and might require adjustements if applied to a different example

% auxiliary files needed: subtightplot.m

clear;
%suptightplot is called (must be included in the same directory)
%for nicer plots
opt = {0.03, 0.03, 0.03} ;
subplot = @(m,n,p) subtightplot(m,n,p,opt{:}); 
figure();

%resol=number of transportation rays and resolution of the integration grid
resol=15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOT 1 OF 2 - not-symmetric case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,1);

%%% NUMBER OF TRANSPORTATION RAYS (the more the rays the more accurate
%%% the final estimate the more time consuming the code)
% upper bound for the domain of a Gaussian distribution
% everything will be restricted to [-M,M]^2
M=4;
%%% INPUTS: COVARIANCE MATRICES of the two Guassian measures %%%
%mean vectors are assumed to be null
mux = [0 0];
muy = [0 0];
global SX SY sx sy;
%upperbound for unbounded numerical integrals 
MMM=6;
%covariance matrix of X,Y. They must be invertible.
SX = [1 0.8; 0.8 2];
SY = [3 1.7; 1.7 1];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOT 2 OF 2 - symmetric case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,2);

%%% NUMBER OF TRANSPORTATION RAYS (the more the rays the more accurate
%%% the final estimate the more time consuming the code)
% upper bound for the domain of a Gaussian distribution
% everything will be restricted to [-M,M]^2

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
hold off;