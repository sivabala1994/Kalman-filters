clc, clear all, close all;
%constants initialised to one.
M = 1;
m = 1;
b = 1;
I = 1;
g = 9.8;
l = 1;

p = I*(M+m)+M*m*l^2; %denominator for the A and B matrices

A = [0      1              0           0;
     0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
     0      0              0           1;
     0 -(m*l*b)/p       m*g*l*(M+m)/p  0];
B = [     0;
     (I+m*l^2)/p;
          0;
        m*l/p];
C = [1 0 0 0;
     0 0 1 0];
D = [0;
     0];

 dt=0.1;
 tfin = 5;
 tspan = 0:dt:tfin;
 
 sysc = ss(A,B,C,D);
 sysd = c2d(sysc,dt);
 
 Ad = sysd.A;
Bd = sysd.B;
Cd = sysd.C;
Dd = sysd.D;

%Simulation
mu = [1 0 0.1 0]';
sigma = eye(4);

u =0;
w1 = normrnd(0,0.1,100,1);
w3 = normrnd(0,0.2,100,1);
w2 = normrnd(0,0.1,100,1);
w4 = normrnd(0,0.3,100,1);
%disturbances
wn= [w1,w2,w3,w4];
%disturbance covariance matrix
Q = [0.1 0 0 0; 
    0 0.2 0 0;
    0 0 0.1 0;
    0 0 0 0.3]; 

%gaussian noise
 v1 = normrnd(0,.2,100,1);
 v3 = normrnd(0,.2,100,1);
 %poisson noise
v1=poissrnd(0.1,100,1);
 v3=poissrnd(0.1,100,1);
vn= [v1,v3]; % as sensor readings are taken ony from two output states
  R = [0.1 0;
      0 0.1];
%  R=cov(vn);
G=eye(4);

H = [1 0;0 1];

%non-linear system simulation
xreal(:,1) = mu;
for i=1 :100
     t= xreal(3,i);
  x= xreal(1,i);
 v= xreal(2,i);
  w= xreal(4,i);
    xreal(:,i+1)=[v;
(4*(0+w^2*sin(t)*0.5 - 0.5*x + cos(t)*sin(t)*0.25))/(4-cos(t)^2);
w;
(4*((-0.5*sin(t)) - 0 - w^2*cos(t)*sin(t)*0.25 + x*cos(t)*0.25))/(4-cos(t)*sin(t))];
    yreal(:,i) = C*(xreal(:,i));
end

%noisy system
xn(:,1) = mu;
for i=1 :100
     xn(:,i+1) = xreal(:,i+1)+wn(i,:)';
     y(:,i) = C*xn(:,i) + H*vn(i,:)';
end



% %% Extended kalman filter


x_p=zeros(101,4);
x_p(1,:)=mu;
s_p=zeros(4,4,101);
s_p(:,:,1)=eye(4);
s_u=zeros(4,4,100);
x_u=zeros(100,4);

for i=1:100
  c(:,:,i)=[x_p(i,1) 0 0 0;
      0 0 x_p(i,3) 0];
  
 s_u(:,:,i)=s_p(:,:,i)-s_p(:,:,i)*c(:,:,i)'*inv(c(:,:,i)*s_p(:,:,i)*c(:,:,i)'+R)*c(:,:,i)*s_p(:,:,i);
 x_u(i,:)=x_p(i,:)+(s_p(:,:,i)*c(:,:,i)'*inv(c(:,:,i)*s_p(i)*c(:,:,i)'+R)*(y(:,i)-(c(:,:,i)*x_p(i,:)')))';

 
 t= x_u(i,3);
  x= x_u(i,1);
 v= x_u(i,2);
  w= x_u(i,4);
 
 a = [                               0, 1,                                                                                                                                                                                                                   0,                                         0;
                2/(cos(t)^2 - 4), 0,                                                                               - (2*w^2*cos(t) + cos(t)^2 - sin(t)^2)/(cos(t)^2 - 4) - (2*cos(t)*sin(t)*(2*sin(t)*w^2 + 0 - 2*x + cos(t)*sin(t)))/(cos(t)^2 - 4)^2,              -(4*w*sin(t))/(cos(t)^2 - 4);
                               0, 0,                                                                                                                                                                                                                   0,                                         1;
 -cos(t)/(4*(cos(t)*sin(t) - 4)), 0, (2*cos(t) + (w^2*cos(t)^2)/4 - (w^2*sin(t)^2)/4 - 0 + (x*sin(t))/4)/(cos(t)*sin(t) - 4) - ((cos(t)^2 - sin(t)^2)*((cos(t)*sin(t)*w^2)/4 + 2*sin(t) + 0 - (x*cos(t))/4))/(cos(t)*sin(t) - 4)^2, (w*cos(t)*sin(t))/(2*(cos(t)*sin(t) - 4))];
 
 
 
 
 

 x_p(i+1,:)=[v;
(4*(0+w^2*sin(t)*0.5 - 0.5*x + cos(t)*sin(t)*0.25))/(4-cos(t)^2);
w;
(4*((-0.5*sin(t)) - 0 - w^2*cos(t)*sin(t)*0.25 + x*cos(t)*0.25))/(4-cos(t)*sin(t))]';


   s_p(:,:,i+1)=a*s_u(:,:,i)*a' + Q;
   yk(:,i)=C*x_u(i,:)';
end


plot(1:100,yreal(2,:),1:100,yk(2,:))
legend('true measurement','EKF estimated measurement'), axis([0 100 -0.5 0.5]),xlabel('time'),ylabel('Position observations')
