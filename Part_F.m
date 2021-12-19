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
c2 = [0 0 1 0 0 0; 0 0 0 0 1 0]; %% theta1, theta2
c3 = [1 0 0 0 0 0; 0 0 0 0 1 0]; %% x(t), theta2
c4 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0]; %% x(t), theta1, theta2
Obs1 = rank([c1' A'*c1' ((A')^2)*c1' ((A')^3)*c1' ((A')^4)*c1' ((A')^5)*c1'])
Obs2 = rank([c2' A'*c2' ((A')^2)*c2' ((A')^3)*c2' ((A')^4)*c2' ((A')^5)*c2'])
Obs3 = rank([c3' A'*c3' ((A')^2)*c3' ((A')^3)*c3' ((A')^4)*c3' ((A')^5)*c3'])
Obs4 = rank([c4' A'*c4' ((A')^2)*c4' ((A')^3)*c4' ((A')^4)*c4' ((A')^5)*c4'])

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

% For C1
L1 = place(A', c1', observer_poles)';
A_obs1 = [(A-B*K_lqr) B*K_lqr; 
        zeros(size(A)) (A-L1*c1)];
B_obs1 = [B;zeros(size(B))];
C_obs1 = [c1 zeros(size(c1))];
ss_obs1 = ss(A_obs1, B_obs1, C_obs1, D)

figure(1)
initial(ss_obs1, initial_state)
figure(2)
step(ss_obs1)

% For C3
L3 = place(A', c3', observer_poles)';
A_obs3 = [(A-B*K_lqr) B*K_lqr; 
        zeros(size(A)) (A-L3*c3)];
B_obs3 = [B;zeros(size(B))];
C_obs3 = [c3 zeros(size(c3))];
ss_obs3 = ss(A_obs3, B_obs3, C_obs3, D)

figure(3)
initial(ss_obs3, initial_state)
figure(4)
step(ss_obs3)

% For C4
L4 = place(A', c4', observer_poles)';
A_obs4 = [(A-B*K_lqr) B*K_lqr; 
        zeros(size(A)) (A-L4*c4)];
B_obs4 = [B;zeros(size(B))];
C_obs4 = [c4 zeros(size(c4))];
ss_obs4 = ss(A_obs4, B_obs4, C_obs4, D)

figure(5)
initial(ss_obs4, initial_state)
figure(6)
step(ss_obs4) 
