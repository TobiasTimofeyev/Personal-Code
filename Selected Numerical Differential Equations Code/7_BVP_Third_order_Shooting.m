% Third order BVP shooting method

clear all;

h=0.02;
u0=[1 0 0]';
v0=[0 1 0]';
w0=[0 0 1]';
xspan=[1,2];  xpnts=xspan(1):h:xspan(2);
%matrix functions for all 3rd-order matrix systems
fu= @(x,Y) [0 1 0; 0 0 1; (1/x^3) (-1/x^2) 0]*Y+[0 0 (log(x)-3)/x^3]'; 
fv= @(x,Y) [0 1 0; 0 0 1; (1/x^3) (-1/x^2) 0]*Y; 
fw= @(x,Y) [0 1 0; 0 0 1; (1/x^3) (-1/x^2) 0]*Y;    
    

%% 2nd order DE Modified Euler method

[xu,Yu] = Tobias_ModEuler_Matrix( fu, u0, xspan, h);

[xv,Yv] = Tobias_ModEuler_Matrix( fv, v0, xspan, h); 

[xw,Yw] = Tobias_ModEuler_Matrix( fw, w0, xspan, h); 

IC=[1/2 1/4]';
A=[Yv(2,length(xpnts)) Yw(2,length(xpnts)); Yv(3,length(xpnts)) Yw(3,length(xpnts))];
theta=A\(IC-Yu(2:3,length(xpnts)));     %solving for appropriate theta values to shoot with

Z=Yu+theta(1)*Yv+theta(2)*Yw;
                                
%% Plotting Results

figure(70201);

% Plot the modified Euler solution:  ---------->
plot(xpnts,Z);
xlabel('x', 'fontsize', 12);   ylabel('Z', 'fontsize', 12);
title('modified Euler shooting','fontsize',14)








                                










