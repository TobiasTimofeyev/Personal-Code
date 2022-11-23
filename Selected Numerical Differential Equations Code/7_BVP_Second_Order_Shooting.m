clear all;

y2=1;             %second BVP boundary value to shoot for
ParamIC=[-10,-11];  %initial theta guesses for finding parameter value
h=0.01;
v0=1;
xspan=[0,2];  xpnts=xspan(1):h:xspan(2);
%matrix functions for all 3rd-order matrix systems
f1v= @(x,y,v) v;
f2v= @(x,y,v) y^2/(2+x); 
    

%% shooting for y2 with 2nd order DE Modified Euler method

tol=true; n=2;
Theta(1)=ParamIC(1);             %initial theta guesses for finding zeros
Theta(2)=ParamIC(2);
%corresponding initial values of F
[~,shoot,~]=   Tobias_ModEuler_2ndOrder( f1v, f2v, v0, Theta(1), xspan, h);
F(1)= (shoot(length(xpnts)))-1;
[~,shoot,~]=   Tobias_ModEuler_2ndOrder( f1v, f2v, v0, Theta(2), xspan, h);
F(2)= (shoot(length(xpnts)))-1;


while tol==true                  %iterating secant method
    
    [~,shoot,~]=   Tobias_ModEuler_2ndOrder( f1v, f2v, v0, Theta(n), xspan, h);
    F(n)= (shoot(length(xpnts)))-y2;
   
    Theta(n+1)=Theta(n)- F(n)/( (F(n)-F(n-1))/(Theta(n)-Theta(n-1)) );
    
    if abs(F((n)))<0.001
        tol=false;
    end
    n=n+1
%     if n>10000      %stop if taking too long 
%         return
%     end
end

Z=shoot;
                 
%% Plotting Results

figure(70301);

% Plot the modified Euler solution:  ---------->
plot(xpnts,Z);
xlabel('x', 'fontsize', 12);   ylabel('Z', 'fontsize', 12);
title('modified Euler shooting','fontsize',14)








                                










