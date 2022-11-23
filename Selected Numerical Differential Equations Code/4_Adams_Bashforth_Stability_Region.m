clear all;

[zR,zI]=meshgrid(-2.1:0.001:0.5,-2.1:0.001:2);        %creates grid of points
z= abs( 0.5*( ( 1 + zR+i*zI + 0.75*(zR+i*zI).^2 ) + sqrt( ( 1 + zR+i*zI + 0.75*(zR+i*zI).^2 ).^2 - ( zR + i*zI ).^2 ) ));  
                                %defining rho_1 function on point matrices
z1= abs( 0.5*( ( 1 + zR+i*zI + 0.75*(zR+i*zI).^2 ) - sqrt( ( 1 + zR+i*zI + 0.75*(zR+i*zI).^2 ).^2 - ( zR + i*zI ).^2 ) ));
                                %defining rho_2 function on point matrices
v=[1,1];                        %specified level == 1
figure(040501);
contour(zR,zI,z,v);             %making the contour
xlabel('h*lambda_R');ylabel('h*lambda_I');
hold on;
contour(zR,zI,z1,v);