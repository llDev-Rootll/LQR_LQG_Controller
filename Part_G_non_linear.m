clear all;
clc;
close all;

%% Defining variables
syms m1 g m2 M l1 l2

Load_Params;

%% Substituting the values of constants

M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 9.81;

A = double(subs(A));
B = double(subs(B));
tspan = 0:0.1:100;
%Enter initial conditions
x0 = 0.5; theta1_0 = deg2rad(10); theta2_0 = deg2rad(10);
initial_state = [ x0, 0, theta1_0, 0, theta2_0, 0];

c1 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];

observer_poles = [-1; -1.5; -2; -2.5; -3; -3.5];
Q = zeros(6);
Q(1,1) = 3000;
Q(2,2) = 0;
Q(3,3) = 1000000;
Q(4,4) = 0;
Q(5,5) = 1000;
Q(6,6) = 0;
R=0.001;
D=0;

K_lqr = lqr(A, B, Q, R);

Vd = 0.02 * eye(size(A));  % Process disturbance
Vn = 0.01; % measurement noise

K_kf = lqr(A', c1', Vd, Vn)';

%% Non-linear Model LQG Response
[t,q1] = ode45(@(t,q)Obs1(t,q,-K_lqr*q,K_kf),tspan,initial_state);
figure();
hold on
plot(t,q1(:,1))
ylabel('state variable')
xlabel('time (sec)')
title('Non-Linear System LQG for condition 1')
legend('x')
hold off

function dQ = Obs1(t,y,F,L)
m1 = 100; m2 = 100; M = 1000; l1 = 20; l2 = 10; g = 9.81;
x = y(1);
dx = y(2);
t1 = y(3);
dt1 = y(4);
t2 = y(5);
dt2 = y(6);
dQ=zeros(6,1);
y1 = [x; 0; 0];
c1 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
sum = L*(y1-c1*y);
dQ(1) = dx + sum(1);
dQ(2) = (F-((m1*sin(t1)*cos(t1))+(m2*sin(t2)*cos(t2)))*g - (l1*m1*(dQ(3)^2)*sin(t1)) - (l2*m2*(dQ(5)^2)*sin(t2)))/(m1+m2+M-(m1*(cos(t1)^2))-(m2*(cos(t2)^2)))+sum(2);
dQ(3) = dt1+sum(3);
dQ(4) = ((cos(t1)*dQ(2)-g*sin(t1))/l1) + sum(4);
dQ(5) = dt2 + sum(5);
dQ(6) = (cos(t2)*dQ(2)-g*sin(t2))/l2 + sum(6);
end