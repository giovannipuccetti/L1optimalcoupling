% Extra material for the paper
% General construction and classes of explicit L1-optimal couplings
% by G. Puccetti ane L. Rüschendorf

% This code computes the trasportation cost given by the coupling described
% in the paper in the case A in Table 6.1,
% and the corresponding dual bound, showing optimality of the coupling construction

% The parameters of the following code are calibrated to the specific case
% Code might not properly work and might require adjustements if applied to a different example

% auxiliary files needed: subtightplot.m / 

clear all
%parameters of the two GAUSSIAN distributions
%null mean
mux = [0 0];
muy = [0 0];
sx=-0.4;
sy=-sx;
%covariance matrix of X,Y
SX = [1 sx; sx 1];
SY = [1 sy; sy 1];
%restriction to [0,M]^2
M=4;

%%% %%% %%% PLOT OF THE DOMAIN OF THE REDUCED DENSITIES %%% %%% %%% %%% 
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
%plot the two diagonals
hold on
plot(u,u,'k','LineWidth',2)
hold on
plot(u,-u,'k','LineWidth',2)
xticks([]);
yticks([]);
axis square
hold on;

%TRANSPORTATION COST by the coupling described in Section 4.4
fun2a = @(x,y) 2*abs(y).*(1/(2*pi*sqrt(1-sx^2))).*(exp((-1/(2*(1-sx^2))).*(x.^2+y.^2+2.*sy.*x.*y))-exp((-1/(2*(1-sx^2))).*(x.^2+y.^2-2.*sy.*x.*y)));
ymin = @(x) -x;
int2a = integral2(fun2a,0,Inf,ymin,0);
format long;
primalbound=4*int2a

%DUAL bound as in (4.15) 
dualbound=2/sqrt(pi)*(sqrt(1+sy)-sqrt(1+sx))