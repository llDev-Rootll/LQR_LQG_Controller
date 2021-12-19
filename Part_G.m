clear all;
clc;
close all;

%% Defining variables
syms m1 g m2 M l1 l2

%% check observability for difference outputs
A = [0 1 0 0 0 0; 0 0 -m1*g/M 0 -m2*g/M 0; 0 0 0 1 0 0; 0 0 -((M*g)+(m1*g))/(M*l1) 0 -g*m2/(M*l1) 0; 0 0 0 0 0 1; 0 0 -m1*g/(M*l2) 0 -((M*g)+(m2*g))/(M*l2) 0];
B = [0; 1/M; 0; 1/(l1*M); 0; 1/(l2*M)]

%% Substituting the values of constants
%%
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 9.81;

A = double(subs(A));
B = double(subs(B));


c1 = [1 0 0 0 0 0]; %% x(t)

observer_poles = [-1; -1.5; -2; -2.5; -3; -3.5];
x0 = 0.5; theta1_0 = deg2rad(10); theta2_0 = deg2rad(10);
initial_state = [ x0, 0, theta1_0, 0, theta2_0, 0, 0, 0, 0,0 ,0 ,0 ];
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

Vd = 0.02*eye(size(A));  % Process disturbance
Vn = 0.01; % measurement noise

K_kf = lqr(A', c1', Vd, Vn)';

ss_lqg = ss([(A-B*K_lqr) B*K_lqr; zeros(size(A)) (A-K_kf*c1)], [B;zeros(size(B))],[c1 zeros(size(c1))], D);
figure(1)
initial(ss_lqg, initial_state)
figure(2)
step(ss_lqg)
