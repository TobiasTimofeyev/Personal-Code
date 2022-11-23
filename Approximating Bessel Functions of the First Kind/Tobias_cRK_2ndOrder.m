

function [x,Y,V]=Tobias_cRK_2ndOrder(f,y0,v0,xspan,h)

x=xspan(1):h:xspan(2); %creates h-steps on xspan interval
Y(1)=y0;                %initial y-value
V(1)=v0;                %initial v-value


for n=1:length(x)-1            %iterating steps
    k1b=h*V(n);                     %my b is y
%     k1a=h*f(Y(n));                %my a is v
    k1a=h*feval(f,x(n),Y(n),V(n));                %my a is v
    
    k2b=h*(V(n)+0.5*k1a);
%     k2a=h*f(Y(n)+0.5*k1b);
    k2a=h*feval(f,x(n)+0.5*h,Y(n)+0.5*k1b,V(n)+0.5*k1a);
    
    k3b=h*(V(n)+0.5*k2a);
%     k3a=h*f(Y(n)+0.5*k2b);
    k3a=h*feval(f,x(n)+0.5*h,Y(n)+0.5*k2b,V(n)+0.5*k2a);
    
    k4b=h*(V(n)+k3a);
%     k4a=h*f(Y(n)+k3b);
    k4a=h*feval(f,x(n)+h,Y(n)+k3b,V(n)+k3a);
    
    Y(n+1)=Y(n)+(1/6)*(k1b+2*k2b+2*k3b+k4b);          %cRK method Y
    V(n+1)=V(n)+(1/6)*(k1a+2*k2a+2*k3a+k4a);          %cRK method V

end

end