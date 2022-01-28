% Extra material for the paper
% General construction and classes of explicit L1-optimal couplings
% by G. Puccetti ane L. Rüschendorf

% This code approximates the trasportation cost given by the coupling described
% in the paper in the case C in Table 6.1 (boundary of reduced measure
% domains is flat)

% The parameters of the following code are calibrated to the specific case
% Code might not properly work and might require adjustements if applied to a different example

% auxiliary files needed: condition3.m /costC.m

%result
% with M=4 and resol=200 bound=0.565354533754756

% with M=4 and resol=300 bound=

clear;

%%% NUMBER OF TRANSPORTATION RAYS (the more the rays the more accurate
%%% the final estimate the more time consuming the code)

%resol=number of transportation rays and resolution of the integration grid
resol=250;
% upper bound for the domain of a Gaussian distribution
% everything will be restricted to [-M,M]^2
M=4;
%%% INPUTS: COVARIANCE MATRICES of the two Guassian measures %%%
%mean vectors are assumed to be null
mux = [0 0];
muy = [0 0];
global xdiag SX SY sx sy d1 d2 s0;
%upperbound for numerical integrals 
MMM=10;
%covariance matrix of X,Y. They must be invertible.
SX = [1.8 -0.4; -0.4 1.8];
SY = [1 0.4; 0.4 1];
dx=det(SX);
dy=det(SY);
%correlation parameters to be used in the density formulas
sx=SX(1,2)/sqrt(SX(1,1)*SX(2,2));
sy=SY(1,2)/sqrt(SY(1,1)*SY(2,2));


%%% %%% %%% %%% TRANSPORTATION RAYS AND  COST %%% %%% %%% %%% %%% %%% %%% 

%rays starting points on the main diagonal
xt=linspace(0,M,resol);
%auxiliary vectors/parameters
dir1=linspace(0,0,resol);
dir2=linspace(0,0,resol);
upp=linspace(0,0,resol);
ddd=linspace(0,0,resol);

s0b=linspace(0,0,resol);
%calculation of boundary points
for t=1:resol
dir1(t)=sqrt(2)/2;
dir2(t)=-sqrt(2)/2;
%point on the upper diagonal
xdiag=xt(t);
%recall directions
d1=sqrt(2)/2;
d2=-sqrt(2)/2;
%calculation of boundary point
s0b(t)=fzero(@condition3,[0,MMM]);
%plot(xdiag+s0b(t)*d1,xdiag+s0b(t)*d2,'ko');
end
s0=max(s0b);

N=resol;
x = linspace(0,M,N);
y = linspace(0,s0,N);
%Alternative method of building the trapz grid for integral calculation
for i=1:N
for j=1:N
    F(i,j)=costC(x(i),y(j));
end
end
format long
4*trapz(y,trapz(x,F,2))

