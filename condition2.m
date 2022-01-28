function s=condition2(rho)
%auxiliary function that returns the difference between the P and Q masses located
%left to R_t, whent R_t is a ray which start from xdiag with angle rho
%used to check condition (6.5) to find optimal direction of transportation
global xdiag;
%global sx sy;
global SX SY sx sy;
%choice of directions
k1=cos(rho);
k2=sin(rho);
%densities of the two measures
f1 = @(t,u) (1/(2*pi*sqrt(SX(1,1)*SX(2,2)*(1-sx^2)))).*(exp((-1/(2*(1-sx^2))).*((t+u.*k1).^2/SX(1,1)+(t+u.*k2).^2/SX(2,2)-2.*sx.*((t+u.*k1)/sqrt(SX(1,1))).*((t+u.*k2)/sqrt(SX(2,2))))));
f2 = @(t,u) (1/(2*pi*sqrt(SY(1,1)*SY(2,2)*(1-sy^2)))).*(exp((-1/(2*(1-sy^2))).*((t+u.*k1).^2/SY(1,1)+(t+u.*k2).^2/SY(2,2)-2.*sy.*((t+u.*k1)/sqrt(SY(1,1))).*((t+u.*k2)/sqrt(SY(2,2))))));
if (rho>5/4*pi)&&(rho<7/4*pi)
ymax = @(t) -2*t/(k1+k2);  
%first part 
intA1=0;
intA2=0;
%second part
intB1 = integral2(f1,0,xdiag,0,ymax);
intB2 = integral2(f2,0,xdiag,0,ymax);
else
%first part    
ymin = @(t) -2*t/(k1+k2); 
%MMM is the upperboudn for unbounded integrals
MMM=20;
intA1 = integral2(f1,-10,0,ymin,Inf);
intA2 = integral2(f2,-10,0,ymin,Inf);
%second part
intB1 = integral2(f1,0,xdiag,0,Inf);
intB2 = integral2(f2,0,xdiag,0,Inf);
end
s=abs(intA1)+abs(intB1)-abs(intA2)-abs(intB2);