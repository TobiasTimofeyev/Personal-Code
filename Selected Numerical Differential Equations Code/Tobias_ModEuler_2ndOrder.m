

function [x,Y,V]=Tobias_ModEuler_2ndOrder(f,y0,v0,xspan,h)

x=xspan(1):h:xspan(2); %creates h-steps on xspan interval
Y(1)=y0;                %initial y-value
V(1)=v0;                %initial v-value

for n=1:length(x)-1            %iterating steps
    Ybar=Y(n)+h*V(n);
    Vbar=V(n)+h*f(Y(n));
    
    V(n+1)=V(n)+(h/2)*(f(Y(n))+f(Ybar));
    
    Y(n+1)=Y(n)+(h/2)*(V(n)+Vbar);       %Modified euler method
end

end