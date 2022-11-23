clear all;

h=0.01;
% h=input('enter a step size h=');     
gamma=0.015;    w0=1;
y0=0; v0=1;
tspan=[0,50];
f= @(x,y,v) -2*gamma*v-w0^2*y; tpnts=tspan(1):h:tspan(2);

w=sqrt(w0^2-gamma^2);
fExact= @(t) (1/w)*exp(-gamma*t).*sin(w*t);
fvExact= @(t) (1/w)*exp(-gamma*t).*(-gamma*sin(w*t)+w*cos(w*t));

EnergyExact= @(t) 0.5*(fvExact(t).^2+w^2*fExact(t).^2);

%% 2nd order Simple Euler

[x1,y1,v1] = Tobias_Euler_2ndOrder( f, y0, v0, tspan, h);   %Calling simple Euler function 
                                %y1 is vector of numerically approximated values
                                
 E1=0.5*(v1.^2+w^2*y1.^2);
% Ham1er=arrayfun(HExact,xpnts1)-Ham1;

%% 2nd order Symplectic Euler (method 2)

[x2,y2,v2] = Tobias_SympEuler_2ndOrder( f, y0, v0, tspan, h);   %Calling simple Euler function 
                                %y1 is vector of numerically approximated values
 E2=0.5*(v2.^2+w^2*y2.^2);
Err1=E1(round(tpnts/h)+1)-EnergyExact(tpnts);
Err2=E2(round(tpnts/h)+1)-EnergyExact(tpnts);
%% plotting

                                
figure(50801)
subplot(2,2,1); plot(y1,v1);
hold on; xlabel('y-values'); ylabel('x-values'); title('Simple Euler');
plot(fExact(tpnts),fvExact(tpnts),':')
axis([ -1.1*max(abs(y1)) 1.1*max(abs(y1)) ...
       -1.1*max(abs(v1)) 1.1*max(abs(v1)) ])
hold off;

subplot(2,2,2); plot(y2,v2);
hold on; xlabel('y-values'); ylabel('x-values'); title('Symplectic Euler 1');
plot(fExact(tpnts),fvExact(tpnts),':')
axis([ -1.1*max(abs(y2)) 1.1*max(abs(y2)) ...
       -1.1*max(abs(v2)) 1.1*max(abs(v2)) ])
hold off;

subplot(2,2,3); plot(tpnts,Err1);
hold on; xlabel('t-values'); ylabel('Error'); title('Energy Error');
plot(tpnts,Err2,'--')
axis([min(tpnts) max(tpnts)... 
    -0.2*max(abs(Err1)) 1.1*max(abs(Err1)) ])
hold off;