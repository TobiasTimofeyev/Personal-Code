%Crank nicholson scheme on HE with homogenous BCs
clear all;
% format long;

h=0.1; 
%need to run these in order
k=0.1; kk=1;
% k=0.05; kk=2;
% k=0.025; kk=3;
% k=0.01; kk=4;
r=k/h^2;

xspan=[0,1]; tspan=[0,0.5];     %setting up grid
xpnts=xspan(1):h:xspan(2); tpnts=tspan(1):k:tspan(2);
M=length(xpnts)-1;

%creating matrix A
A=diag( ones(1,M-2) ,-1) -  2*diag( ones(1,M-1) )  + diag( ones(1,M-2) ,1);

u0= sin(pi*xpnts);      %initial condition

B1=eye(M-1)+(r/2)*A;  %matrices needed for step calculations
B2=eye(M-1)-(r/2)*A;

U=zeros(length(xpnts), length(tpnts)); %making matrix for computed values
U(:,1)=u0;

for n=1:length(tpnts)-1         %time iterations of CN
    rvec=B1 *U(2:end-1,n);
    U(2:end-1,n+1)= thomas( diag(B2,-1), diag(B2), diag(B2,1), rvec);
    %note: leaving boundaries as zero due to BCs.
end

%line colors and styles
line_color(1,:)=[1 0 0];   % red
line_color(2,:)=[0 0 0];   % black 
line_color(3,:)=[0 0 1];   % blue
line_color(4,:)=[0 0.8 0.2];   % greenish
line_style=char('-','--','-.',':');

%plotting solution
figure(130401)
hold on;
title('Plot at final time step');
plot(xpnts,U(:,end), 'color',line_color(kk,:),'linestyle',line_style(kk,:),'linewidth',2)
xlabel('x'); ylabel('t'); zlabel('u');
legend('k=0.1', 'k=0.5', 'k=0.025', 'k=0.01','Location', 'Best')
hold off;

%plotting error
figure(130402)
UEndExact= sin(pi*xpnts) *exp(-pi^2*0.5);
hold on;
title('error at final time step');
plot(xpnts,UEndExact'-U(:,end), 'color',line_color(kk,:),'linestyle',line_style(kk,:),'linewidth',1)
xlabel('x'); ylabel('t'); zlabel('error');
legend('k=0.1', 'k=0.5', 'k=0.025', 'k=0.01','Location', 'Best')
hold off;


%point (0.5,0.5) error
disp('the error at (t,x)=(0.5,0.5) is ')
Point_Error= sin(pi*0.5)*exp(-pi^2*0.5)-U( round(0.5/h)+1 , round(0.5/k)+1 )

%3D plot
figure(130403)
[X,T]=meshgrid(xpnts,tpnts);
surf(X,T,U')
xlabel('x'); ylabel('t'); zlabel('u');





