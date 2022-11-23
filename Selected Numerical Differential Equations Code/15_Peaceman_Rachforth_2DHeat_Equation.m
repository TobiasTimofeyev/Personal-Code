%simplicit explicit method for 2D heat equation with nonhomogeneous BCs
clear all;

%setting up values and steps
ymax=2.5;
xspan=[0,1]; yspan=[0,ymax]; tspan=[0,0.5];
h=0.1; k=0.002; r=k/h^2;
xpnts=xspan(1):h:xspan(2); ypnts=yspan(1):h:yspan(2); tpnts=tspan(1):k:tspan(2);
[X,Y]=meshgrid(xpnts,ypnts);
M=length(xpnts)-1; L=length(ypnts)-1;

U0=@(x,y) 10*sin(pi*x)*sin(pi*y/ymax) + sin(2*pi*y/ymax)*(1-x) + cos(2*pi*y/ymax)*x;
uexact=@(x,y) 10*sin(pi*x)*sin(pi*y/ymax)*exp(-(1+1/ymax^2)*pi^2*tspan(2)) +...
    (sin(2*pi*y/ymax)*(1-x)+cos(2*pi*y/ymax)*x)*exp(-(2*pi/ymax)^2*tspan(2));

%initial condition and exact solution.
%Assuming that the IC satisfies BCs
Uold=zeros(length(ypnts),length(xpnts)); Unew=zeros(length(ypnts),length(xpnts)); 
for l=1:length(ypnts)
    for m=1:length(xpnts)
        Uold(l,m)=U0(xpnts(m),ypnts(l));
        UE(l,m)=uexact(xpnts(m),ypnts(l));
    end
end

%time iteration
for p=1:length(tpnts)-1
    %BCs
    Unew(1,:)= xpnts * exp(-(2*pi/ymax)^2*tpnts(p+1));   
    Unew(length(ypnts),:)= xpnts * exp(-(2*pi/ymax)^2*tpnts(p+1));  
    Unew(:,1)= sin(2*pi*ypnts/ymax)*exp(-(2*pi/ymax)^2*tpnts(p+1));
    Unew(:,length(xpnts))= cos(2*pi*ypnts/ymax)*exp(-(2*pi/ymax)^2*tpnts(p+1));
    
    for l=2:length(ypnts)-1
        for m=2:length(xpnts)-1
            U=Uold;
            Unew(l,m)=U(l,m)+r*(U(l,m+1)-2*U(l,m)+U(l,m-1))+r*(U(l+1,m)-2*U(l,m)+U(l-1,m));
        end
    end
%     surf(X,Y,Unew)
%     pause;
    Uold=Unew;
end

figure(150201)
surf(X,Y,Unew)
title('numerical Solution')
xlabel('x'); ylabel('y'); zlabel('z');

figure(150202)
surf(X,Y,UE-Unew)
title('Error')
xlabel('x'); ylabel('y'); zlabel('z');