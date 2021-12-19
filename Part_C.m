clear all;
clc;
close all;

%% Defining variables
syms m1 g m2 M l1 l2

%% check observability for difference outputs
A = [0 1 0 0 0 0; 0 0 -m1*g/M 0 -m2*g/M 0; 0 0 0 1 0 0; 0 0 -((M*g)+(m1*g))/(M*l1) 0 -g*m2/(M*l1) 0; 0 0 0 0 0 1; 0 0 -m1*g/(M*l2) 0 -((M*g)+(m2*g))/(M*l2) 0];
B = [0; 1/M; 0; 1/(l1*M); 0; 1/(l2*M)]
ctrb = [B A*B ((A)^2)*B ((A)^3)*B ((A)^4)*B ((A)^5)*B];
det(ctrb)