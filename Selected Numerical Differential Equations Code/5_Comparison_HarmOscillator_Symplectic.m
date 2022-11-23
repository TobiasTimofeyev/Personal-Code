clear all;

h1=0.2; h2=0.5;
y0=0; v0=1;
xspan=[0,1000];  xpnts1=xspan(1):h1:xspan(2);  xpnts2=xspan(1):h2:xspan(2);
f= @(y) -y;

fExact= @(x) sin(x);
fvExact= @(x) cos(x);
HExact= @(x) 0.5;

%% 2nd order DE Verlet1

[x1,y1,v1] = Tobias_Verlet1_2ndOrder( f, y0, v0, xspan, h1);   %Calling simple Euler function 
                                %y1 is vector of numerically approximated values
                                
Ham1=0.5*v1.^2+0.5*y1.^2;
Ham1er=arrayfun(HExact,xpnts1)-Ham1;


%% 2nd order DE Modified Euler method

[x2,y2,v2] = Tobias_ModEuler_2ndOrder( f, y0, v0, xspan, h1);   %Calling cRK function 
                                %y2 is vector of numerically approximated values

Ham2=0.5*v2.^2+0.5*y2.^2;
Ham2er=arrayfun(HExact,xpnts1)-Ham2;
                                

%% 2nd order DE cRK method

[x3,y3,v3] = Tobias_cRK_2ndOrder( f, y0, v0, xspan, h2);   %Calling this implicit solution function 
                                %y3 is vector of numerically approximated values

Ham3=0.5*v3.^2+0.5*y3.^2;
Ham3er=arrayfun(HExact,xpnts2)-Ham3;
                                
%% 2nd order DE central difference method

[x4,y4,v4] = Tobias_Center_2ndOrder( f, y0, v0, xspan, h1);   %Calling this implicit solution function 
                                %y4 is vector of numerically approximated values

                                
Ham4=0.5*v4.^2+0.5*y4.^2;
Ham4er=arrayfun(HExact,xpnts1)-Ham4;

%% 2nd order DE ode45

z=@(t,y) vectorF(t,y);
options=odeset('AbsTol',0.002);
[x5_A,y5_A] = ode45( z, xspan, [y0,v0], options);   %Calling this implicit solution function 
                                %y5_A is 2d vector of numerically approximated values
                                
options=odeset('AbsTol',0.003);
[x5_B,y5_B] = ode45( z, xspan, [y0,v0], options);   %Calling this implicit solution function 
                                %y5_A is 2d vector of numerically approximated values

Ham5_A=0.5*(y5_A(:,2).^2+(y5_A(:,1)).^2);
Ham5er_A=arrayfun(HExact,x5_A)-Ham5_A;

Ham5_B=0.5*(y5_B(:,2).^2+(y5_B(:,1)).^2);
Ham5er_B=arrayfun(HExact,x5_B)-Ham5_B;
                                
                                
                                
                                
%% Plotting Results

figure(50601);
% Plot the Verlet solution:  ---------->
subplot(3,2,1);  plot(y1, v1);
axis([ -1.1*max(abs(y1)) 1.1*max(abs(y1)) ...
       -1.1*max(abs(v1)) 1.1*max(abs(v1)) ])
axis('equal')
xlabel('y', 'fontsize', 12);   ylabel('v', 'fontsize', 12);
title('Verlet-1','fontsize',14)
% Plot the modified Euler solution:  ---------->
subplot(3,2,2);  plot(y2, v2);
axis([ -1.1*max(abs(y2)) 1.1*max(abs(y2)) ...
       -1.1*max(abs(v2)) 1.1*max(abs(v2)) ])
axis('equal')
xlabel('y', 'fontsize', 12);   ylabel('v', 'fontsize', 12);
title('modified Euler','fontsize',14)
% Plot the cRK solution:  ---------->
subplot(3,2,3);  plot(y3, v3);
axis([ -1.1*max(abs(y3)) 1.1*max(abs(y3)) ...
       -1.1*max(abs(v3)) 1.1*max(abs(v3)) ])
axis('equal')
xlabel('y', 'fontsize', 12);   ylabel('v', 'fontsize', 12);
title('classical Runge--Kutta','fontsize',14)
% Plot the central-difference solution:  ---------->
subplot(3,2,4);  plot(y4, v4);
axis([ -1.1*max(abs(y4)) 1.1*max(abs(y4)) ...
       -1.1*max(abs(v4)) 1.1*max(abs(v4)) ])
axis('equal')
xlabel('y', 'fontsize', 12);   ylabel('v', 'fontsize', 12);
title('central-difference','fontsize',14)
% Plot the ode45 solution with AbsTol = 0.002:  ---------->
subplot(3,2,5);  plot(y5_A(:,1), y5_A(:,2));
axis([ -1.01 1.01 -1.01 1.01 ])
axis('equal')
xlabel('y', 'fontsize', 12);   ylabel('v', 'fontsize', 12);
title('ode45 with Tol=0.002','fontsize',14)
%Plot the ode45 solution with AbsTol = 0.003:  ---------->
subplot(3,2,6);  plot(y5_B(:,1), y5_B(:,2));
axis([ -1.01 1.01 -1.01 1.01 ])
axis('equal')
xlabel('y', 'fontsize', 12);   ylabel('v', 'fontsize', 12);
title('ode45 with Tol=0.003','fontsize',14)


% figure(343)
% 
% plot(xpnts1,Ham1er)


figure(50602);
% Plot the Verlet Hamiltonian error:  ---------->
subplot(3,2,1);  plot(xpnts1, Ham1er); axis([ 0 1.1*max(abs(xpnts1)) -1.1*max(abs(Ham1er)) 1.1*max(abs(Ham1er)) ])
xlabel('x', 'fontsize', 12);   ylabel('H-Error', 'fontsize', 12);
title('Verlet-1-H','fontsize',14)
% Plot the modified Euler Hamiltonian error:  ---------->
subplot(3,2,2);  plot(xpnts1, Ham2er); axis([ 0 1.1*max(abs(xpnts1)) -1.1*max(abs(Ham2er)) 1.1*max(abs(Ham2er)) ])
xlabel('x', 'fontsize', 12);   ylabel('H-Error', 'fontsize', 12);
title('modified Euler-H','fontsize',14)
% Plot the cRK Hamiltonian error:  ---------->
subplot(3,2,3);  plot(xpnts2, Ham3er); axis([ 0 1.1*max(abs(xpnts1)) -1.1*max(abs(Ham3er)) 1.1*max(abs(Ham3er)) ])
xlabel('x', 'fontsize', 12);   ylabel('H-Error', 'fontsize', 12);
title('classical Runge--Kutta-H','fontsize',14)
% Plot the central-difference Hamiltonian error:  ---------->
subplot(3,2,4);  plot(xpnts1, Ham4er); axis([ 0 1.1*max(abs(xpnts1)) -1.1*max(abs(Ham1er)) 1.1*max(abs(Ham1er)) ])
xlabel('x', 'fontsize', 12);   ylabel('H-Error', 'fontsize', 12);
title('central-difference-H','fontsize',14)
% Plot the ode45 Hamiltonian error with AbsTol = 0.002:  ---------->
subplot(3,2,5);  plot(x5_A, Ham5er_A); axis([ 0 1.1*max(abs(x5_A))  -1.1*max(abs(Ham5er_A)) 1.1*max(abs(Ham5er_A)) ])
xlabel('x', 'fontsize', 12);   ylabel('H-Error', 'fontsize', 12);
title('ode45 Tol=0.002-H','fontsize',14)
%Plot the ode45 Hamiltonian error with AbsTol = 0.003:  ---------->
subplot(3,2,6);  plot(x5_B, Ham5er_B); axis([ 0 1.1*max(abs(x5_B)) -1.1*max(abs(Ham5er_B)) 1.1*max(abs(Ham5er_B)) ])
xlabel('x', 'fontsize', 12);   ylabel('H-Error', 'fontsize', 12);
title('ode45 Tol=0.003-H','fontsize',14)











                                










