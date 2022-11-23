

function [x,Y,V]=Tobias_Center_2ndOrder(f,y0,v0,xspan,h)

x=xspan(1):h:xspan(2); %creates h-steps on xspan interval
Y(1)=y0;                %initial y-values
Y(2)=y0+h*v0+(h^2/2)*f(y0);                

V(1)=v0;                %initial v-value

for n=2:length(x)-1            %iterating steps
    Y(n+1)=h^2*f(Y(n))+2*Y(n)-Y(n-1); %obtaining y-values with center difference
end

for n=2:length(x)-2
    V(n+1)=(Y(n+2)-Y(n))/(2*h);
end
V(length(x))=(Y(length(x))-Y(length(x)-1))/(h) - (h/2)*(Y(length(x)));

end