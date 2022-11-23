%Finite element method collocation method
clear all;

M=40;   %number of points/ equations
xspan=[0,1]; h=(xspan(2)-xspan(1))/(M+1); xpnts=xspan(1):h:xspan(2);

%solving for y after variable transform to z=y-x
%these are the terms in DE z''+Q(x)*z=R(x)
R=@(x) (2*x-4)./((1+x).^2);     
Q=@(x) (-2)./((1+x).^2);
phi=@(j,x) sin(j.*pi.*x);
ddphi=@(j,x) -(j^2)*(pi^2)*sin(j*pi*x);

%Exact solution
yExact=@(x) 2*x./(1+x); 
A=zeros(M,M);

for i=2:(M+1)
    r(i-1)=R(xpnts(i));
    for j=2:(M+1)
       A(i-1,j-1)= feval(ddphi,j-1,xpnts(i)) + Q(xpnts(i))*feval(phi,j-1,xpnts(i)); 
    end
end
r=r';

% xpnts=xspan(1):0.01:xspan(2);
%FEM raw solution and then final solution after sum and variable change
c=A\r;
for m=1:length(xpnts)
    basis=feval(phi,1:M, xpnts(m))';
    z(m)= sum(c.*basis);
end
yvals=z+xpnts;

yE= yExact(xpnts);   %exact solution points
figure(90101)
hold on
plot(xpnts,yvals,'b'); plot(xpnts,yE,'r');
xlabel('x'); ylabel('y');
hold off


figure(90102)
hold on
plot(xpnts,yE-yvals)
xlabel('x'); ylabel('error');
hold off
