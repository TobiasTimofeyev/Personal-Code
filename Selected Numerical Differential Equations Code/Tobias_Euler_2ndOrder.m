
function [x,Y,V]=Tobias_Euler_2ndOrder(f,y0,v0,xspan,h)

x=xspan(1):h:xspan(2); %creates h-steps on xspan interval
Y(1)=y0;                %initial y-value
V(1)=v0;                %initial v-value

for n=1:length(x)-1            %iterating steps
    V(n+1)=V(n)+h*(feval(f,x(n),Y(n),V(n)));
    
    Y(n+1)=Y(n)+h*(V(n));       %Simple euler method
end

end