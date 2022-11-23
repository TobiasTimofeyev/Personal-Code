clear all;

[zR,zI]=meshgrid(-3:0.1:1,-2:0.1:2);        %creates grid of points
z= (1+zR+0.5*(zR.^2-zI.^2)).^2+(zI+zR.*zI).^2;  
                                %defining function on point matrices
v=[1,1];                        %specified level == 1
figure(040101);
contour(zR,zI,z,v);             %making the contour
xlabel('h*lambda_R');ylabel('h*lambda_I');