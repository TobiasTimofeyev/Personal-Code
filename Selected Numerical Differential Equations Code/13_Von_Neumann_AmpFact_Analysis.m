%question, plotting amplification factors
%When doing von neuman analysis on a numerical Scheme, an amplification
%factor, rho, is found. This determines the growth of error in the problem
%This code was used in the analysis of this amplification factor  for
%different parameter values in the problem.
clear all;

h=0.05;

betah=0:0.01:pi;

theta1=1; theta2=0.5;

r1=1/(2*h); r2=1/(4*h);

%both r-values for theta=1
for n=1:length(betah)
 rho1Theta1(n)= (1-(1-theta1)*4*r1*sin(betah(n)/2)^2) / (1+theta1*4*r1*sin(betah(n)/2)^2);
end
 
 for n=1:length(betah)
 rho2Theta1(n)= (1-(1-theta1)*4*r2*sin(betah(n)/2)^2) / (1+theta1*4*r2*sin(betah(n)/2)^2);
 end
 
 %both r-values for theta=1/2
 for n=1:length(betah)
 rho1Theta2(n)= (1-(1-theta2)*4*r1*sin(betah(n)/2)^2) / (1+theta2*4*r1*sin(betah(n)/2)^2);
 end
 
 for n=1:length(betah)
 rho2Theta2(n)= (1-(1-theta2)*4*r2*sin(betah(n)/2)^2) / (1+theta2*4*r2*sin(betah(n)/2)^2);
 end
 %exact r1
 rhoTheta1E=exp(-r1*(betah).^2);
 %exact r2
 rhoTheta2E=exp(-r2*(betah).^2);
 
 %line colors and styles
line_color(1,:)=[1 0 0];   % red
line_color(2,:)=[0 0 0];   % black 
line_color(3,:)=[0 0 1];   % blue
line_color(4,:)=[0 0.8 0.2];   % greenish
line_color(5,:)=[0 0.8 0.2];   % greenish
line_color(6,:)=[0 0.8 0.2];   % greenish
line_style=char('-','--','-.',':');
 
 figure(1305011)
 hold on
 plot(betah,rho1Theta1, 'color',[1 0 0],'linestyle', '--','linewidth',1)
 plot(betah,rho2Theta1, 'color',[0 0 0],'linestyle', '-.','linewidth',1)
 plot(betah,rhoTheta1E, 'color',[0 0 1],'linestyle', '-','linewidth',1)
 plot(betah,rho1Theta2, 'color',[1 1 0],'linestyle', '--','linewidth',1)
 plot(betah,rho2Theta2, 'color',[0 1 0],'linestyle', '-.','linewidth',1)
 plot(betah,rhoTheta2E, 'color',[0 1 1],'linestyle', '-','linewidth',1)
 legend('theta=1, r=1/2h', 'theta=1, r=1/4h', 'r=1/2h Exact', 'theta=1/2, r=1/2h','theta=1/2, r=1/4h', 'r=1/4h Exact')
 xlabel('beta h'); ylabel('rho');
 hold off

 
 
 