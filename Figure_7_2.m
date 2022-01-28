% Extra material for the paper
% General construction and classes of explicit L1-optimal couplings
% by G. Puccetti ane L. Rüschendorf

% This code produces Figure 7.2  which illustrates Directions of optimal transport (left) 
% and optimal dual function (right) between discrete versions of the two Gaussian marginals 
% as in case A in Table 6.1. 
% Discrete Gaussian distributions are obtained via N = 500 simulations. 
% On the left picture we plot the domains of the reduced densities, whereas on the right the exact dual function (7.3) 
% is superimposed on the discretized dual one.

% auxiliary files needed: subtightplot.m

% requires (easy) installation of CVX solver for Matlab; see http://cvxr.com/cvx/

clear all
%suptightplot is called (the .m file must be included in the same directory)
%for nicer plots
opt = {0.03, 0.03, 0.03} ;
subplot = @(m,n,p) subtightplot(m,n,p,opt{:}); 
figure();
rng(202);
%N is the cardinality of marginal measures
N=500;
%marginal probabilites are simple uniform on N points
prob=ones(1,N)/N;
%GENERATION OF MARGINALS (set of points for bivariate normal)
%null mean
mux = [0 0];
muy = [0 0];
%correlation parameters
sx=-0.40;
sy=0.40;
%covariance matrix of X,Y
SX = [1 sx; sx 1];
SY = [1 sy; sy 1]; 
%Gaussian distributions limited to [-M,M]^2
M=4;
% we fix the random seed for reproducibility
rng('default')
%N simulations from the marginals
X = mvnrnd(mux,SX,N);
Y = mvnrnd(muy,SY,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LEFT PLOT - discretized coupling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,1);
%values grid
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
plot(30:30, 'k.-','MarkerSize',5), axis([-M M -M M]);
hold on
image([-M,M],[-M,M],fun)
colormap gray
axis square
hold on;

%%%% P R I M A L P R O B L E M %%%%%

cmat=zeros(N,N);
        for i=1:N
            for j=1:N
                %here L1 distance but can be changed to different cost
                cmat(i,j)=norm(X(i,:)-Y(j,:));
            end
        end
%LP begins subject to CVX notation
cvx_begin
    variable x(N,N);
    OBJ=sum(sum(cmat.*x));
    tic;
    minimize(OBJ)
    subject to
    %marginal constraints
    sum(x,1)==prob
    sum(x,2)==transpose(prob)
    %non negativity constraints
        x>=0
%LPend
cvx_end;
%computational time for one instance
comp=toc;

%figure representing the transportation between points by red line
[A,B] = find(x>0.00001);
for i=1:N 
plot([X(A(i),1), Y(B(i),1)], [X(A(i),2), Y(B(i),2)],'k','LineWidth',2)
axis([-M M -M M])
axis square
hold on
end

%plot the two diagonals to illustrate the corresponding sectors

x = linspace(-4,4,100);
y = x;
plot(x,x,'--k')
hold on
plot(x,-x,'--k')
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RIGHT PLOT - discretized dual function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,2);

%%%% D U A L P R O B L E M %%%%%
%cost vector CMAT
CMAT=zeros(N^2,1);
h=1;
        for i=1:N
            for j=1:N
                CMAT(h)=norm(X(i,:)-Y(j,:));
                h=h+1;
            end
        end
%constraint matrix CNST
e=ones(1,N);
CNST=[kron(eye(N),e)' repmat(eye(N),N,1)];
%LP begins
cvx_begin
%first N variables refer to first dimension; second N variables to second
%dimensions
    variable x(2*N);
    OBJ=sum(sum(x))/N;
    tic;
    maximize(OBJ)
    subject to
    %dual constraints
    CNST*x<=CMAT
%LPend
cvx_end;
%computational time for one instance
comp=toc;

%now we look at the dual variables.
f=x(1:N);
g=x((N+1):(2*N));
%dual functions can be given up to a constant by the LP
f=f-mean(f);
g=g-mean(g);

%plot real dual funztion in 3d
x = linspace(-4,4,25);
y = x';
z = (abs(x.^1-y.^1)-abs(x.^1+y.^1))/2;
s=surf(x,y,z)
colormap('gray');
axis([-4 4 -4 4 -4 4])
alpha(s,.3)
hold on;
%3d plot of the discretized dual functions
scatter3(X(:,1),X(:,2),f,[], 'Black','LineWidth',4)
hold on
scatter3(Y(:,1),Y(:,2),-g, [], 'Black','LineWidth',4)
hold off;