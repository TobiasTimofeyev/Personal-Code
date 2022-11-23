%2D Heat equation solved using peaceman-Rachforth
clear all;

%setting up values and steps
ymax=2.5;
xspan=[0,1]; yspan=[0,ymax]; tspan=[0,0.5];
h=0.1; k=h/4; r=k/h^2;
xpnts=xspan(1):h:xspan(2); ypnts=yspan(1):h:yspan(2); tpnts=tspan(1):k:tspan(2);
[X,Y]=meshgrid(xpnts,ypnts);
M=length(xpnts)-1; L=length(ypnts)-1;

U0=@(x,y) 10*sin(pi*x)*sin(pi*y/ymax) + sin(2*pi*y/ymax)*(1-x) + cos(2*pi*y/ymax)*x;
uexact=@(x,y) 10*sin(pi*x)*sin(pi*y/ymax)*exp(-(1+1/ymax^2)*pi^2*tspan(2)) +...
    (sin(2*pi*y/ymax)*(1-x)+cos(2*pi*y/ymax)*x)*exp(-(2*pi/ymax)^2*tspan(2));

%initial condition and exact solution.
%Assuming that the IC satisfies BCs
Uold=zeros(length(ypnts),length(xpnts));
for l=1:length(ypnts)
    for m=1:length(xpnts)
        Uold(l,m)=U0(xpnts(m),ypnts(l));
        UE(l,m)=uexact(xpnts(m),ypnts(l));
    end
end

Unew=zeros(length(ypnts),length(xpnts)); 
bstarvec=zeros(M-1,1); bvec=zeros(L-1,1);

%time iteration
for p=1:length(tpnts)-1
    %BCs for solution at next time step
    Unew(1,:)= xpnts * exp(-(2*pi/ymax)^2*tpnts(p+1));   
    Unew(length(ypnts),:)= xpnts * exp(-(2*pi/ymax)^2*tpnts(p+1));  
    Unew(:,1)= sin(2*pi*ypnts/ymax)*exp(-(2*pi/ymax)^2*tpnts(p+1));
    Unew(:,length(xpnts))= cos(2*pi*ypnts/ymax)*exp(-(2*pi/ymax)^2*tpnts(p+1));
    
    %Vectors determining Ustar at boundaries
    Gvec0=0.5*(Uold(2:L ,1) + Unew(2:L ,1)) + (r/4)*((Uold(3:L+1,1) -2*Uold(2:L,1) + Uold(1:L-1,1))...
                 - (Unew(3:L+1,1) -2*Unew(2:L,1) + Unew(1:L-1,1)));
    Gvec1=0.5*(Uold(2:L ,M+1) + Unew(2:L ,M+1)) + (r/4)*((Uold(3:L+1,M+1) -2*Uold(2:L,M+1) + Uold(1:L-1,M+1))...
                 - (Unew(3:L+1,M+1) -2*Unew(2:L,M+1) + Unew(1:L-1,M+1)));
    
    %implementing BCs to Ustar
    Ustar=zeros(length(ypnts),length(xpnts));
    Ustar(2:L,1)=Gvec0; Ustar(2:L,M+1)=Gvec1;
    
    for l=2:L
        bstarvec(1)=Gvec0(l-1); bstarvec(M-1)=Gvec1(l-1);   %making BC dependent vector
        cvec= Uold(l, 2:M)' + (r/2)*( Uold(l+1, 2:M)' - 2*Uold(l, 2:M)'+ Uold(l-1, 2:M)' ) + (r/2)*bstarvec;
        Ustar(l, 2:M)=thomas( -r/2*ones(M-2,1) , (1+r)*ones(M-1,1) , -r/2*ones(M-2,1) , cvec );
    end
    
    
    for m=2:M
        bvec(1)=Unew(1,m); bvec(L-1)=Unew(L+1,m);   %making BC dependent vector
        cvec=Ustar(2:L, m) + (r/2)*( Ustar(2:L, m+1) - 2*Ustar(2:L, m)+ Ustar(2:L, m-1) ) +(r/2)*bvec;
        Unew(2:L, m)=thomas( -r/2*ones(L-2,1) , (1+r)*ones(L-1,1) , -r/2*ones(L-2,1) , cvec );
    end
%     figure(2)
%     surf(X,Y,Ustar)
%     figure(1)
%     surf(X,Y,Unew)
%     pause;
    Uold=Unew;
end

figure(150801)
surf(X,Y,Unew)
title('Numerical Solution')
xlabel('x'); ylabel('y'); zlabel('z');

figure(150802)
surf(X,Y,UE-Unew)
title('Error')
xlabel('x'); ylabel('y'); zlabel('z');