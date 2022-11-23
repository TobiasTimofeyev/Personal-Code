%Crank nicholson scheme on HE with variable coeffients
clear all;

h=0.02; 
%need to run these in order
k=h/8;
r=k/h^2;

%coefficients for equation we're solving: u_tt=(alpha*u_x)_x+beta*u
alpha=@(x) 1+x; 
beta=@(x,t) (-x^2/(1+t^2));

xspan=[0,1]; tspan=[0,0.5];     %setting up grid
xpnts=xspan(1):h:xspan(2); tpnts=tspan(1):k:tspan(2);
M=length(xpnts)-1;

u0= sin(pi*xpnts);      %initial condition

U=zeros(length(xpnts), length(tpnts)-3); %making matrix for computed values
U(:,1)=u0;


%defining a base matrix for A and B to come from (no beta term on main diagonal)
for m=1:M-2
   a(m)= -0.5*r* 0.5*( alpha(xpnts(m+1)) + alpha(xpnts(m+2)) ); %subdiagonal
   
   b(m)= 1+ 0.5*r* 0.5*( alpha(xpnts(m+1)) + alpha(xpnts(m+2)) ) + ...
      0.5*r* 0.5*( alpha(xpnts(m)) + alpha(xpnts(m+1)) ); %main diagonal
  
   c(m)= -0.5*r* 0.5*( alpha(xpnts(m+1)) + alpha(xpnts(m+2)) ); %superdiagonal
end
b(M-1)= 1+ 0.5*r*0.5*( alpha(xpnts(M)) + alpha(xpnts(M+1)) ) + ...
       0.5*r*0.5*( alpha(xpnts(M-1)) + alpha(xpnts(M)) ); %last diagonal term
% c=[0 c]; a=[a 0];   %appending zeros to accomodate spdiags syntax
a=a'; b=b'; c=c';

D=diag(a,-1) + diag(b,0) + diag(c,1);
   
% D=spdiags(a,-1,M-1,M-1) + spdiags(b,0,M-1,M-1) + spdiags(c,1,M-1,M-1);
    

    

for n=1:length(tpnts)-2         %time iterations of CN
    
    avec=0;bvec=0;%reset
    for m=1:length(xpnts)-2
        avec(m)= -0.5*k*beta(xpnts(m+1),tpnts(n+2));
        bvec(m)= 0.5*k*beta(xpnts(m+1),tpnts(n+1));
    end
    avec=avec'; bvec=bvec';
    
%     A= D + spdiags(avec,0,M-1,M-1);
    A= D + diag(avec,0);
    
%     B= -D + 2*spdiags(ones(M-1,1),0,M-1,M-1) + spdiags(bvec,0,M-1,M-1);
    B= -D + 2*eye(M-1) + diag(bvec,0);
    
    
    rvec=B *U(2:end-1,n);
    
    U(2:end-1,n+1)= A\rvec;
    
%     figure(140201)
%     plot(xpnts,U(:,n+1), 'color','r','linewidth',2)
%     pause;
%     
    %note: leaving right boundary as zero due to BCs.
end

% %plotting solution
figure(140201)
hold on;
title('Plot at final time step (t=0.5)');
plot(xpnts,U(:,end), 'color','r','linewidth',2)
xlabel('x'); ylabel('u');
hold off;

% 3D plot
figure(140202)
[X,T]=meshgrid(xpnts,tpnts(1:end-1));
surf(X,T,U')
xlabel('x'); ylabel('t'); zlabel('u');

