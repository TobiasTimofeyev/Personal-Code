% BVP Shooting method
% Solving BVP with IVPs
clear all;

h=0.1;
u0=1; uprime0=0;
v0=0; vprime0=1;
xspan=[0,2];  xpnts=xspan(1):h:xspan(2);
f= @(x,y,z) z;              %u and v share first equation in first order system
g1= @(x,y,z) 3*x-x*z+3*y;   %second equation in first order system for u
g2= @(x,y,z) -x*z+3*y;      %second equation in first order system for homog. v

%% 2nd order DE Modified Euler method

[xu,yu,zu] = Tobias_ModEuler_2ndOrder( f, g1, u0, uprime0, xspan, h);

[xv,yv,zv] = Tobias_ModEuler_2ndOrder( f, g2, v0, vprime0, xspan, h); 

theta=(5-yu(length(xpnts)))/yv(length(xpnts));
w=yu+theta*yv;
                                
%% Plotting Results

figure(70101);

% Plot the modified Euler solution:  ---------->
plot(xpnts,w);
xlabel('x', 'fontsize', 12);   ylabel('y', 'fontsize', 12);
title('modified Euler shooting','fontsize',14)








                                










