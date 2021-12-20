clear all;
clc;
close all;

%% Symbols
syms F M m1 m2 l1 l2 g;
syms x xd xdd theta1 theta1d theta1dd theta2 theta2d thetadd;

%% Defining state variables
q = [x xd theta1 theta1d theta2 theta2d];
u = [F];

%% Defining the eqns - Non-Linear system

% double derivative of x
xdd = (F - m1*g*cos(theta1)*sin(theta1) -  m2*g*cos(theta2)*sin(theta2) - m1*l1*(theta1d^2)*sin(theta1) - m2*l2*(theta2d^2)*sin(theta2))/(M + m1*(sin(theta1))^2 + m2*(sin(theta2))^2);

% double derivative of theta1
theta1dd = (xdd*cos(theta1) - g*sin(theta1))/l1;

% double derivative of theta2
theta2dd = (xdd*cos(theta2) - g*sin(theta2))/l2;

%% Non-linear Model
Model = [xd xdd theta1d theta2dd theta2d theta2dd];

%% Linearization
%% 
% Jacobian
J = jacobian(Model, q);

q_e = [0 0 0 0 0 0];


A = subs(J, q, q_e);
B = [0; 1/M; 0; 1/(l1*M); 0; 1/(l2*M)];
C = eye(size(A));
D=0;

C_m = [B A*B (A^2)*B (A^3)*B (A^4)*B (A^5)*B];
r = rank(C_m);

M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 9.81;
A = double(subs(A))
B = double(subs(B));
eigs(A)
Q = zeros(6);
Q(1,1) = 3000;
Q(2,2) = 0;
Q(3,3) = 1000000;
Q(4,4) = 0;
Q(5,5) = 1000;
Q(6,6) = 0;
R=0.001;

K_lqr = lqr(A, B, Q, R);

Ac = A - B*K_lqr;
eigs(Ac)

open_sys = ss(A, B, C, D);
x0 = 0.5; theta1_0 = deg2rad(10); theta2_0 = deg2rad(10);

initial_state = [ x0, 0, theta1_0, 0, theta2_0, 0];

state_feedback_ss = ss(Ac, B, C, D);

initial( state_feedback_ss, initial_state);