function z=condition4(s)
%auxiliary function
global xdiag;
%global sx sy;
global SX SY sx sy d1 d2 int2 int1 s1 s0 delta1 delta2;
%conditional densities
f1 = @(u) (abs(d2-d1+u*(delta1*d2-delta2*d1)).*(1/(2*pi*sqrt(SX(1,1)*SX(2,2)*(1-sx^2)))).*(exp((-1/(2*(1-sx^2))).*((xdiag+u.*d1).^2/SX(1,1)+(xdiag+u.*d2).^2/SX(2,2)-2.*sx.*((xdiag+u.*d1)/sqrt(SX(1,1))).*((xdiag+u.*d2)/sqrt(SX(2,2)))))));
f2 = @(u) (abs(d2-d1+u*(delta1*d2-delta2*d1)).*(1/(2*pi*sqrt(SY(1,1)*SY(2,2)*(1-sy^2)))).*(exp((-1/(2*(1-sy^2))).*((xdiag+u.*d1).^2/SY(1,1)+(xdiag+u.*d2).^2/SY(2,2)-2.*sy.*((xdiag+u.*d1)/sqrt(SY(1,1))).*((xdiag+u.*d2)/sqrt(SY(2,2)))))));
%conditional densities for the reduced measures
g1 = @(u) max((f1(u)-f2(u)),0);
g2 = @(u) max((f2(u)-f1(u)),0);
z = integral(g2,0,s1)/int2-integral(g1,s0,s)/int1;
