
%semi-implicit solution for nonlinear pde 
clear all;

h=0.1; 
%need to run these in order
k=h/4;
r=k/h^2;
gamma=0.25;
b0=0; b1=2; %fixed BCs

xspan=[-5,5]; tspan=[0,3];     %setting up grid
xpnts=xspan(1):h:xspan(2); tpnts=tspan(1):k:tspan(2);
M=length(xpnts)-1;

a=@(x) x-1;

u0= (4/pi)*atan(exp(7*xpnts));      %initial condition

U=zeros(length(xpnts), length(tpnts)-1); %making matrix for computed values
U(:,1)=u0;

%% Initial Step using explicit method

for m=2:length(xpnts)-1
    U(m,2)=U(m,1)+a(U(m,1))*r*h*0.5*(U(m+1,1)-U(m-1,1))+gamma*r*(U(m+1,1)-2*U(m,1)+U(m-1,1));
end
U(end,2)=b1;    %end nonzero 

%% Main iteration

%main diagonals don't change in time
aMain=(1+gamma*r)*ones(M-1,1);
bMain=(1-gamma*r)*ones(M-1,1);

BC=zeros(M-1,1); %setting up BC vector

for n=2:length(tpnts)-1
    Uhalf=1.5*U(:,n)-0.5*U(:,n-1);  %U^(n+1/2) approximation
    
    %setting up upper/lower diags for matrices (time dependent)
    for m=1:M-2
        aUpper(m)=-gamma*(r/2)-a(Uhalf(m+1))*(r*h/4);
        aLower(m)=-gamma*(r/2)+a(Uhalf(m+2))*(r*h/4);
        
        bUpper(m)=gamma*(r/2)+a(Uhalf(m+1))*(r*h/4);
        bLower(m)=gamma*(r/2)-a(Uhalf(m+2))*(r*h/4);
    end
    
    B=diag(bUpper,1) + diag(bMain,0) + diag(bLower,-1);
    
    BC(end)=b1* (gamma*r+a(Uhalf(end-1))*0.5*r*h);    %adjusting time-dependent BC vector
    
    %Solving linear system for next time row
    bvec=B*U(2:end-1,n) + BC;
    
    U(2:end-1,n+1)=thomas(aLower,aMain,aUpper,bvec);
    U(end,n+1)=b1;  %adding BC
end


uE= 1+tanh(xpnts/(2*gamma));              %exact solution


%plotting Numerical and Exact solution
figure(140301)
hold on;
title('Plot at final time step (t=3)');
plot(xpnts,U(:,end), 'color','r','linewidth',2)
plot(xpnts,uE, 'color','g','linewidth',2,'linestyle',':')
xlabel('x'); ylabel('u');
hold off;

%plotting Error
figure(140302)
hold on;
title('Error at final time step (t=3)');
plot(xpnts,uE-U(:,end)', 'color','r','linewidth',2)
xlabel('x'); ylabel('error');
hold off;


% % 3D plot
 figure(140303)
 [X,T]=meshgrid(xpnts,tpnts);
 surf(X,T,U')
 xlabel('x'); ylabel('t'); zlabel('u')
 



