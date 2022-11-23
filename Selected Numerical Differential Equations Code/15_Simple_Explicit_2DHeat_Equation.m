%simple explicit method for 2D heat equation
clear all;

%setting up values and steps
ymax=2.5;
xspan=[0,1]; yspan=[0,ymax]; tspan=[0,0.5];
h=0.1; k=0.002; r=k/h^2;
xpnts=xspan(1):h:xspan(2); ypnts=yspan(1):h:yspan(2); tpnts=tspan(1):k:tspan(2);
[X,Y]=meshgrid(xpnts,ypnts);
M=length(xpnts)-1; L=length(ypnts)-1;

uexact=@(x,y) sin(pi*x) * sin(pi*y/ymax) * exp(-(1+ymax^(-2))*pi^2*tspan(2));

%initial condition and exact solution
Uold=ones(length(ypnts),length(xpnts));   
for l=1:length(ypnts)
    for m=1:length(xpnts)
        Uold(l,m)=sin(pi*xpnts(m))*sin((pi*ypnts(l))/ymax);
        UE(l,m)=uexact(xpnts(m),ypnts(l));
    end
end

Unew=zeros(length(ypnts),length(xpnts)); 

%time iteration
for p=1:length(tpnts)-1
    %Boundary Conditions. Redundant in this case. This is how you might
    %otherwise change them
    Unew(1,:)= 0;   
    Unew(length(ypnts),:)= 0;
    Unew(:,1)= 0;
    Unew(:,length(xpnts))= 0;
    
    for l=2:length(ypnts)-1
        for m=2:length(xpnts)-1
            %Iteration of numerical scheme
            U=Uold;
            Unew(l,m)=U(l,m)+r*(U(l,m+1)-2*U(l,m)+U(l,m-1))+r*(U(l+1,m)-2*U(l,m)+U(l-1,m));
        end
    end
%     surf(X,Y,Unew)
%     pause;
    Uold=Unew;
end

figure(150101)
surf(X,Y,Unew)
title('Numerical Solution')
xlabel('x'); ylabel('y'); zlabel('z');

figure(150102)
surf(X,Y,UE-Unew)
title('Error')
xlabel('x'); ylabel('y'); zlabel('z');