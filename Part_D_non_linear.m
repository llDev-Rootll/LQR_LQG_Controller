clear all
clc

syms F M m1 m2 l1 l2 g;
syms x xd xdd theta1 theta1d theta1dd theta2 theta2d thetadd;

q = [x xd theta1 theta1d theta2 theta2d];
u = [F];
xdd_num = F-(g*m1*(sin(2*theta1))/2)-(g*m2*(sin(2*theta2))/2)- (m1*l1*(theta1d*theta1d)*sin(theta1)) - (m2*l2*(theta2d*theta2d)*sin(theta2));
xdd_den = M + m1*sin(theta1*theta1)+m1*sin(theta1*theta1);
xdd = xdd_num/xdd_den;
theta1dd = (xdd*cos(theta1) - g*sin(theta1))/l1;
theta2dd = (xdd*cos(theta2) - g*sin(theta2))/l2;
Model = [xd xdd theta1d theta2dd theta2d theta2dd];

J = jacobian(Model, q);
q_e = [0 0 0 0 0 0];
A = subs(J, q, q_e);
B = [0; 1/M; 0; 1/(l1*M); 0; 1/(l2*M)];
M = 1000; m1 = 100; m2 = 100; l1 = 20; l2 = 10; g = 9.81;
A=[0 1 0 0 0 0; 0 0 -(m1*g)/M 0 -(m2*g)/M 0; 0 0 0 1 0 0; 0 0 -((M+m1)*g)/(M*l1) 0 -(m2*g)/(M*l1) 0; 0 0 0 0 0 1; 0 0 -(m1*g)/(M*l2) 0 -(g*(M+m2))/(M*l2) 0];
B=[0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];
Q = zeros(6);
Q(1,1) = 1000;
Q(2,2) = 0;
Q(3,3) = 1000000;
Q(4,4) = 0;
Q(5,5) = 1000000;
Q(6,6) = 0;
R = 10;
y0 = [5; 0; deg2rad(30); 0; deg2rad(60); 0];
tspan = 0:0.1:500;
[t1,y1] = ode45(@(t,y) odefcn(t,y,A,B,Q,R,g,M,m1,m2,l1,l2), tspan, y0);
plot(t1,y1)
legend x x-dot theta1 theta1-dot theta2 theta2-dot
ylabel('state variables')
xlabel('time (sec)')
grid on

function dydt = odefcn(t,y,A,B,Q,R,g,M,m1,m2,l1,l2)
K_lqr = lqr(A, B, Q, R);
F = -K_lqr*y;
dydt=zeros(6,1);
dydt(1) = y(2);
dydt(2)=(F-(g/2)*(m1*sind(2*y(3))+m2*sind(2*y(5)))-(m1*l1*(y(4)^2)*sind(y(3)))-(m2*l2*(y(6)^2)*sind(y(5))))/(M+m1*((sind(y(3)))^2)+m2*((sind(y(5)))^2));
dydt(3)= y(4);
dydt(4)= (dydt(2)*cosd(y(3))-g*(sind(y(3))))/l1'; 
dydt(5)= y(6); 
dydt(6)= (dydt(2)*cosd(y(5))-g*(sind(y(5))))/l2;
end
