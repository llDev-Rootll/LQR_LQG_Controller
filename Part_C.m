clear all;
clc;
close all;

syms m1 g m2 M l1 l2

Load_Params;
ctrb = [B A*B ((A)^2)*B ((A)^3)*B ((A)^4)*B ((A)^5)*B];
det(ctrb)