%auxiliary function to integrate to compute the cost in case C
function z=costC(t,s)
global SX SY sx sy d1 d2 s0 MMM;
%functions for conditional densities
f1 = @(u) (sqrt(2).*(1/(2*pi*sqrt(SX(1,1)*SX(2,2)*(1-sx^2)))).*(exp((-1/(2*(1-sx^2))).*((t+u.*d1).^2/SX(1,1)+(t+u.*d2).^2/SX(2,2)-2.*sx.*((t+u.*d1)/sqrt(SX(1,1))).*((t+u.*d2)/sqrt(SX(2,2)))))));
f2 = @(u) (sqrt(2).*(1/(2*pi*sqrt(SY(1,1)*SY(2,2)*(1-sy^2)))).*(exp((-1/(2*(1-sy^2))).*((t+u.*d1).^2/SY(1,1)+(t+u.*d2).^2/SY(2,2)-2.*sy.*((t+u.*d1)/sqrt(SY(1,1))).*((t+u.*d2)/sqrt(SY(2,2)))))));
%conditional densities of the two reduced measures
g1 = @(u) max((f1(u)-f2(u)),0);
g2 = @(u) max((f2(u)-f1(u)),0);
condition5 = @(y) (integral(g2,0,s)-integral(g1,s0,y));
z=abs(fzero(condition5,[s0,MMM])-s)*g2(s);
