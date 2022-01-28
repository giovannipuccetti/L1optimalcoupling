function s=condition2b(rho)
%auxiliary function that returns the difference between the P and Q masses located
%left to R_t, whent R_t is a ray which start from xdiag with angle rho
%used to check condition (6.5) to find optimal direction of transportation%global sx sy;
global SX SY sx sy MMM xell yell;
%choice of directions
k1=cos(rho);
k2=sin(rho);
%point on the diagonal
xdiag=xell+k1*(xell-yell)/(k2-k1);
if xdiag<0
    s=-1;
else
%densities of the two measures
f1 = @(t,u) (1/(2*pi*sqrt(SX(1,1)*SX(2,2)*(1-sx^2)))).*(exp((-1/(2*(1-sx^2))).*((t+u.*k1).^2/SX(1,1)+(t+u.*k2).^2/SX(2,2)-2.*sx.*((t+u.*k1)/sqrt(SX(1,1))).*((t+u.*k2)/sqrt(SX(2,2))))));
f2 = @(t,u) (1/(2*pi*sqrt(SY(1,1)*SY(2,2)*(1-sy^2)))).*(exp((-1/(2*(1-sy^2))).*((t+u.*k1).^2/SY(1,1)+(t+u.*k2).^2/SY(2,2)-2.*sy.*((t+u.*k1)/sqrt(SY(1,1))).*((t+u.*k2)/sqrt(SY(2,2))))));
%second part
intB1 = integral2(f1,xdiag,MMM,0,MMM);
intB2 = integral2(f2,xdiag,MMM,0,MMM);
s=abs(intB1)-abs(intB2);
end