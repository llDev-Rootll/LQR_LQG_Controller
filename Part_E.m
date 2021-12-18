
clear all;
clc;
close all;

%% Defining variables
syms m1 g m2 M L1 L2

%% check observability for difference outputs
A = [0 1 0 0 0 0; 0 0 -m1*g/M 0 -m2*g/M 0; 0 0 0 1 0 0; 0 0 -((M*g)+(m1*g))/(M*L1) 0 -g*m2/(M*L1) 0; 0 0 0 0 0 1; 0 0 -m1*g/(M*L2) 0 -((M*g)+(m2*g))/(M*L2) 0];
B = [0; 1/M; 0; 1/(L1*M); 0; 1/(L2*M)]
c1 = [1 0 0 0 0 0]; %% x(t)
c2 = [0 0 1 0 0 0; 0 0 0 0 1 0]; %% theta1, theta2
c3 = [1 0 0 0 0 0; 0 0 0 0 1 0]; %% x(t), theta2
c4 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0]; %% x(t), theta1, theta2
Obs1 = rank([c1' A'*c1' ((A')^2)*c1' ((A')^3)*c1' ((A')^4)*c1' ((A')^5)*c1'])
Obs2 = rank([c2' A'*c2' ((A')^2)*c2' ((A')^3)*c2' ((A')^4)*c2' ((A')^5)*c2'])
Obs3 = rank([c3' A'*c3' ((A')^2)*c3' ((A')^3)*c3' ((A')^4)*c3' ((A')^5)*c3'])
Obs4 = rank([c4' A'*c4' ((A')^2)*c4' ((A')^3)*c4' ((A')^4)*c4' ((A')^5)*c4'])