% Extra material for the paper
% General construction and classes of explicit L1-optimal couplings
% by G. Puccetti ane L. Rüschendorf

% This code compute the dual bounds 
% for cases A to F in Table 6.1

clear all
%fix random seed
rng(1);
% % %CASE A 
rho=0.4;
%dual sharp bound as in (4.14)
dual_A= 2/sqrt(pi)*((sqrt(1+rho))-(sqrt(1-rho)));

% % % % %CASE B
SX=[1 -0.4; -0.4 1];
SY=[1 0.8; 0.8 1];
%computation the dual bound in (4.15)
kp=sqrt(2/pi)*(sqrt(SX(1,1)+2*SX(1,2)+SX(2,2))-sqrt(SY(1,1)+2*SY(1,2)+SY(2,2)));
km=sqrt(2/pi)*(sqrt(SX(1,1)-2*SX(1,2)+SX(2,2))-sqrt(SY(1,1)-2*SY(1,2)+SY(2,2)));
dual1=sqrt((kp^2+km^2)/2);
a=sqrt(kp^2/2/(kp^2+km^2));
b=sqrt(km^2/2/(kp^2+km^2));

%computation the dual bound in (4.16) 
%computation of the norm of the two vectors by simulations
N=10^8;
X=transpose(mvnrnd([0,0],SX,N));
Y=transpose(mvnrnd([0,0],SY,N));
dual2=abs(sum(vecnorm(X)/N)-sum(vecnorm(Y)/N));
%choose the best dual bound
dual_B=max(dual1,dual2);

% %CASE C
SX=[1.8 -0.4; -0.4 1.8];
SY=[1 0.4; 0.4 1];
%computation the dual bound in (4.15)
kp=sqrt(2/pi)*(sqrt(SX(1,1)+2*SX(1,2)+SX(2,2))-sqrt(SY(1,1)+2*SY(1,2)+SY(2,2)));
km=sqrt(2/pi)*(sqrt(SX(1,1)-2*SX(1,2)+SX(2,2))-sqrt(SY(1,1)-2*SY(1,2)+SY(2,2)));
dual1=sqrt((kp^2+km^2)/2);
a=sqrt(kp^2/2/(kp^2+km^2));
b=sqrt(km^2/2/(kp^2+km^2));

%computation the dual bound in (4.16)  
%computation of the norm of the two vectors by simulations
N=10^8;
X=transpose(mvnrnd([0,0],SX,N));
Y=transpose(mvnrnd([0,0],SY,N));
dual2=abs(sum(vecnorm(X)/N)-sum(vecnorm(Y)/N));
%choose the best dual bound
dual_C=max(dual1,dual2);

% % % %CASE D
SX=[2 -0.4; -0.4 2];
SY=[1 0.4; 0.4 1];
%computation the dual bound in (4.15)
kp=sqrt(2/pi)*(sqrt(SX(1,1)+2*SX(1,2)+SX(2,2))-sqrt(SY(1,1)+2*SY(1,2)+SY(2,2)));
km=sqrt(2/pi)*(sqrt(SX(1,1)-2*SX(1,2)+SX(2,2))-sqrt(SY(1,1)-2*SY(1,2)+SY(2,2)));
dual1=sqrt((kp^2+km^2)/2);
a=sqrt(kp^2/2/(kp^2+km^2));
b=sqrt(km^2/2/(kp^2+km^2));

%computation the dual bound in (4.16)  
%computation of the norm of the two vectors by simulations
N=10^8;
X=transpose(mvnrnd([0,0],SX,N));
Y=transpose(mvnrnd([0,0],SY,N));
dual2=abs(sum(vecnorm(X)/N)-sum(vecnorm(Y)/N));
%choose the best dual bound
dual_D=max(dual1,dual2);

% %CASE E
SX=[1.3 0.4; 0.4 1.3];
SY=[1 0.4; 0.4 1];
%computation the dual bound in (4.15) 
kp=sqrt(2/pi)*(sqrt(SX(1,1)+2*SX(1,2)+SX(2,2))-sqrt(SY(1,1)+2*SY(1,2)+SY(2,2)));
km=sqrt(2/pi)*(sqrt(SX(1,1)-2*SX(1,2)+SX(2,2))-sqrt(SY(1,1)-2*SY(1,2)+SY(2,2)));
dual1=sqrt((kp^2+km^2)/2);
a=sqrt(kp^2/2/(kp^2+km^2));
b=sqrt(km^2/2/(kp^2+km^2));

%computation the dual bound in (4.16) 
%computation of the norm of the two vectors by simulations
N=10^8;
X=transpose(mvnrnd([0,0],SX,N));
Y=transpose(mvnrnd([0,0],SY,N));
dual2=abs(sum(vecnorm(X)/N)-sum(vecnorm(Y)/N));
%choose the best dual bound
dual_E=max(dual1,dual2);

%CASE F 
%for case F dual bound can be computed analytically via (4.6)
dual_F=(sqrt(2)-1)*sqrt(pi/2);

sprintf('%0.4f',dual_A)
sprintf('%0.4f',dual_B)
sprintf('%0.4f',dual_C)
sprintf('%0.4f',dual_D)
sprintf('%0.4f',dual_E)
sprintf('%0.4f',dual_F)