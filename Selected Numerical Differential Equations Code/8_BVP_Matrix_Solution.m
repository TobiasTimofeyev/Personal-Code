%This code demonstrates a numerical solution to a BVP as a matrix equation.
%As each node is dependent on surrounding nodes, this method lets us
%resolve all dependencies at once.
clear all;

h=0.025; %step size
y0=0; beta=2; B1=1; B2=2;
xspan=[0,1];    xpnts=xspan(1):h:xspan(2);
N=round((xspan(2)-xspan(1))/h);      %matrix dimension. 'round' makes it an integer

exact=@(x) 2*x./(1+x);
Yexact= exact(xpnts);

%P,Q,R functions Defining 2nd order BVP Equation 
P=@(x) 0; 
Q=@(x) -2./((1.+x).^2); 
R=@(x) -4./((1.+x).^2); 

%constructing the diagonals of the tridiagonal matrix
a(1:N-2)=1.-(h/2)*P(xpnts(3:N)); 
b(1:N-1)=-1*(2.-h^2*Q(xpnts(2:N)));
c(1:N-1)=1.+(h/2)*P(xpnts(2:N));
a(N-1)=2;
b(N)=-( 2 - h^2*Q(xpnts(N+1)) + 2*h*(B1/B2)*( 1 + (h/2)*P(xpnts(N+1)) ));
%constructing vector r for Ay=r
r(1)= h^2*R(xpnts(2)) - (1-(h/2)*P(xpnts(2)))*y0;
r(2:N-1)=h^2*R(xpnts(3:N));
r(N)=h^2*R(xpnts(N+1)) - 2*h*(beta/B2)*(1+(h/2)*P(xpnts(N+1)));

Y(1)=0;
Y(2:N+1)=thomas(a,b,c,r);

figure(80501)
plot(xpnts,Yexact-Y)
xlabel('x'); ylabel('error');

figure(80502)
plot(xpnts, Y)
xlabel('x'); ylabel('y');