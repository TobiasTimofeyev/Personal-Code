

function [x,Y]=Tobias_ModEuler_Matrix(f,y0,xspan,h)

x=xspan(1):h:xspan(2); %creates h-steps on xspan interval
Y(1:length(y0),1)=y0;                %initial y-value vector

for n=1:length(x)-1            %iterating steps
    Ybar=Y(:,n)+h*feval(f,x(n),Y(:,n));
    
    Y(:,n+1)=Y(:,n)+(h/2)*(feval(f,x(n),Y(:,n))+feval(f,x(n),Ybar));  %matrix modified euler
end

end