
%Numerically approximating bessel function of the first kind
%This was part of a project in which I worked on approximating Bessel
%Functions of the first kind
%The difficulty in doing so lies in the fact that 
clear all;


h=0.01;
xspan=[0,10];  xpnts=xspan(1):h:xspan(2);

%The method of approximation involves approximating the solution of the
%Bessel equation a step away from the singular point through taylor
%expansions, then continuing with classical Runge-Kutta approximation

% 'a' is the eiegenvalue that characterizes a instance of the Bessel equation

%stepping away for a=0
% a=0;
% y0=1; v0=0;
% y1=y0+h*v0;
% v1=v0+h*(-y1/ (2-(a^2)/2));



% stepping away for a=1
% a=1;
% y0=0; v0=0.5;
% y1=y0+h*v0;
% v1=v0+h*(-y1/ (2-(a^2)/2));



% stepping away for a=2
% a=2;
% y0=0; v0=0;
% v0prime=0.25;        %providing z''
% y1=y0+h*v0+ h^2/2 * v0prime; 
% v1=v0+h*v0prime;

 

% stepping away for a=3
a=3;
y0=0; v0=0;


%A key point about the Bessel equation is that constant multiples of
%solutions are solutions. Thus it is up to us to choose a scaling factor
%that matches the conventional solution as closely as possible.
c=1/8; %Choice of z'''

y1=y0+h*v0*+h^2*0+ (h^3/6)*c + h^4*0 + (h^5/120)*(-20/(25-a^2))*c;
v1=v0+h*0+ (h^2/2)*c + h^3*0 + (h^4/24)*(-20/(25-a^2))*c;  



f= @(x,y,v) -(1/x)*v-(1-a^2/x^2)*y; %Bessel Equation



[x3,y3,v3] = Tobias_cRK_2ndOrder( f, y1, v1, [xspan(1)+h,xspan(2)], h);   %Calling this implicit solution function 
%                                 %y3 is vector of numerically approximated values

% [x3,y3,v3] = Tobias_Euler_2ndOrder( f, y1, v1, [xspan(1)+h,xspan(2)], h);   %Calling this implicit solution function 
                                %y3 is vector of numerically approximated values
x3=[0, x3]; y3=[y0, y3];  v3=[v0, v3];                             
                                
%% Plotting Results

% Plot the solution:  ---------->
hold on;
plot(xpnts,y3,'r');
% plot(xpnts,besselj(0,xpnts)-y3);
plot(xpnts,besselj(3,xpnts),'b');
xlabel('x', 'fontsize', 12);   ylabel('Z', 'fontsize', 12);
% title('cRK','fontsize',14)
legend('cRK','Exact')
hold off;

maxError=max(besselj(0,xpnts)-y3)







                                










