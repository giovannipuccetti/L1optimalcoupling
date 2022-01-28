% Extra material for the paper
% General construction and classes of explicit L1-optimal couplings
% by G. Puccetti ane L. Rüschendorf

% This code approximates the trasportation cost given by the coupling described
% in the paper in the case F in Table 6.1 

%In this case the coupling is described in Section 4.3 and
%the transportation cost in analytical:

%covariance matrix of X,Y of the form aI_d and bI_d.
%SX = [a 0; 0 a];
%SY = [b 0; 0 b];
a=1;
b=2;
%optimal transportation cost
sqrt(2)*gamma(3/2)*(sqrt(b)-sqrt(a))