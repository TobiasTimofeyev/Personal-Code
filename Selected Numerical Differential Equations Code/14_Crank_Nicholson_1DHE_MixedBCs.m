%Crank nicholson scheme on HE with mixed BCs
clear all;

h=0.02; 
%need to run these in order
k=h/2;
r=k/h^2;

p=1; q=pi;  %constants

xspan=[0,1]; tspan=[0,2];     %setting up grid
xpnts=xspan(1):h:xspan(2); tpnts=tspan(1):k:tspan(2);
M=length(xpnts)-1;

%creating matrix A
A= (-r/2)*diag( ones(1,M-1) ,-1) +  (1+r)*diag( ones(1,M) )  + (-r/2)*diag( ones(1,M-1) ,1);
A(1,1)=1+r*(1-h*p); A(1,2)=-r;

%creating matrix B
B= (r/2)*diag( ones(1,M-1) ,-1) +  (1-r)*diag( ones(1,M) )  + (r/2)*diag( ones(1,M-1) ,1);
B(1,1)=1-r*(1-h*p); B(1,2)=r;

%creating vector b
bvec=zeros(M,1); bvec(1)=-r*h*2*q;


u0= sin(pi*xpnts);      %initial condition

U=zeros(length(xpnts), length(tpnts)); %making matrix for computed values
U(:,1)=u0;

for n=1:length(tpnts)-1         %time iterations of CN
    rvec=B *U(1:end-1,n) + bvec;
    U(1:end-1,n+1)= thomas( diag(A,-1), diag(A), diag(A,1), rvec);
    %note: leaving right boundary as zero due to BCs.
end

%plotting solution
figure(140101)
hold on;
title('Plot at final time step (t=2)');
plot(xpnts,U(:,end), 'color','r','linewidth',2)
xlabel('x'); ylabel('u');
hold off;

%estimating boundary condition at x=0
LBC=U(1,end) + ( (U(2,end)-U(1,end))/(h) - (h/2)*( U(1,end)-2*U(2,end)+U(3,end) )/h^2 )

% 3D plot
figure(140102)
[X,T]=meshgrid(xpnts,tpnts);
surf(X,T,U')
xlabel('x'); ylabel('t'); zlabel('u');





