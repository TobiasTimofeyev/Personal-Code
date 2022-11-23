%discretized heat equation with Q4 modified method
%This is the graph of numerical scheme for the heat equation where the von
%neuman analysis would tell us that that the error will dominate the
%solution and grow.
clear all;

%setup of steps and stepsizes
h=0.05; k=0.004; r=k/h^2;
xspan=[0,1]; tspan=[0,0.5]; 
xpnts=xspan(1):h:xspan(2); tpnts=tspan(1):k:tspan(2);

%exact solution for reference
yExact=@(x,t) sin(pi*x).*exp(-pi^2*t);
yEF=yExact(xpnts,tspan(2));

%initial condition
ICeq=@(x) sin(pi*x);
U=zeros(length(xpnts),length(tpnts));
U(1:length(xpnts),1)=ICeq(xpnts);
U(1,:)=0; U(end,:)=0;   %boundary conditions

Ubar=zeros(length(xpnts)-2);
%boundary conditions for the supplementary Ubar, assuming BCs are fixed.
%if BCs were time dependent, then these would have to change 
Ubar(1)=U(1,1); Ubar(end)=U(end,1); 

for n=1:length(tpnts)-1  %time steps
    
    for m=2:length(xpnts)-1
       Ubar(m)= r*U(m+1,n) + (1-2*r)*U(m,n) + r*U(m-1,n);
    end
    
    for m=2:length(xpnts)-1  %position steps
         %new method from Q4 that uses Ubar values calculated in previous loop        
          U(m,n+1)= U(m,n) + (r/2)*( (U(m+1,n) - 2*U(m,n) + U(m-1,n)) + (Ubar(m+1) - 2*Ubar(m) + Ubar(m-1)));  
    end
%     pause %plotting at each step
%     hold on;
%     figure(120401)
%     plot(xpnts,U(:,n))
%     xlabel('x'); ylabel('t');
%     hold off;
end


%error plot for t=0.5
figure(120402)
hold on;
plot(xpnts,yEF-U(:,end)');
xlabel('x'); ylabel('error');
hold off;

%3D plot
figure(120403)
[X,T]=meshgrid(xpnts,tpnts);
surf(X,T,U')
xlabel('x'); ylabel('t'); zlabel('u');




