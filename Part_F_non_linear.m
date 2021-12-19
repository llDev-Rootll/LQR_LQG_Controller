clear all

%% Defining variables
syms m1 g m2 M l1 l2
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 9.81;
%q0 = [0.5 0 deg2rad(10) 0 deg2rad(10) 0];
initial_state = [ 0.5, 0, deg2rad(10), 0, deg2rad(10), 0];
tspan = 0:0.1:100;

%% Observability Check
A = [0 1 0 0 0 0; 0 0 -m1*g/M 0 -m2*g/M 0; 0 0 0 1 0 0; 0 0 -((M*g)+(m1*g))/(M*l1) 0 -g*m2/(M*l1) 0; 0 0 0 0 0 1; 0 0 -m1*g/(M*l2) 0 -((M*g)+(m2*g))/(M*l2) 0];
B = [0; 1/M; 0; 1/(l1*M); 0; 1/(l2*M)];
c1 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
c2 = [0 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
c3 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 0];
c4 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
observer_poles = [-1;-1.5;-2;-2.5;-3;-3.5];
Obs1_check = rank([c1' A'*c1' ((A')^2)*c1' ((A')^3)*c1' ((A')^4)*c1' ((A')^5)*c1']);
Obs2_check = rank([c2' A'*c2' ((A')^2)*c2' ((A')^3)*c2' ((A')^4)*c2' ((A')^5)*c2']);
Obs3_check = rank([c3' A'*c3' ((A')^2)*c3' ((A')^3)*c3' ((A')^4)*c3' ((A')^5)*c3']);
Obs4_check = rank([c4' A'*c4' ((A')^2)*c4' ((A')^3)*c4' ((A')^4)*c4' ((A')^5)*c4']);


L1 = place(A', c1', observer_poles)';
L3 = place(A', c3', observer_poles)';
L4 = place(A', c4', observer_poles)';


%%
%% Non-linear Model Observer Response
[t,q1] = ode45(@(t,q)Obs1(t,q,1,L1),tspan,initial_state);
figure();
hold on
plot(t,q1(:,1))
ylabel('state variables')
xlabel('time (sec)')
title('Non-Linear System Observer for condition 1')
legend('x')
hold off

[t,q3] = ode45(@(t,q)Obs3(t,q,1,L3),tspan,initial_state);
figure();
hold on
plot(t,q3(:,1))
plot(t,q3(:,5))
ylabel('state variables')
xlabel('time (sec)')
title('Non-Linear System Observer for condition 3')
legend('x','theta_2')
hold off

[t,q4] = ode45(@(t,q)Obs4(t,q,1,L4),tspan,initial_state);
figure();
hold on
plot(t,q4(:,1))
plot(t,q4(:,3))
plot(t,q4(:,5))
ylabel('state variables')
xlabel('time (sec)')
title('Non-Linear System Observer for condition 4')
legend('x','theta_1','theta_2')
hold off

function dQ = Obs1(t,y,F,L)
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 9.81;
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

function dQ = Obs3(t,y,F,L)
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 9.81;
x = y(1);
dx = y(2);
t1 = y(3);
dt1 = y(4);
t2 = y(5);
dt2 = y(6);
dQ=zeros(6,1);
y3 = [x; 0; t2];
c3 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 0];
sum = L*(y3-c3*y);
dQ(1) = dx + sum(1);
dQ(2) = (F-((m1*sin(t1)*cos(t1))+(m2*sin(t2)*cos(t2)))*g - (l1*m1*(dQ(3)^2)*sin(t1)) - (l2*m2*(dQ(5)^2)*sin(t2)))/(m1+m2+M-(m1*(cos(t1)^2))-(m2*(cos(t2)^2)))+sum(2);
dQ(3) = dt1+sum(3);
dQ(4) = ((cos(t1)*dQ(2)-g*sin(t1))/l1) + sum(4);
dQ(5) = dt2 + sum(5);
dQ(6) = (cos(t2)*dQ(2)-g*sin(t2))/l2 + sum(6);
end

function dQ = Obs4(t,y,F,L)
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 9.81;
x = y(1);
dx = y(2);
t1 = y(3);
dt1 = y(4);
t2 = y(5);
dt2 = y(6);
dQ=zeros(6,1);
y4 = [x; t1; t2];
c4 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
sum = L*(y4-c4*y);
dQ(1) = dx + sum(1);
dQ(2) = (F-((m1*sin(t1)*cos(t1))+(m2*sin(t2)*cos(t2)))*g - (l1*m1*(dQ(3)^2)*sin(t1)) - (l2*m2*(dQ(5)^2)*sin(t2)))/(m1+m2+M-(m1*(cos(t1)^2))-(m2*(cos(t2)^2)))+sum(2);
dQ(3) = dt1+sum(3);
dQ(4) = ((cos(t1)*dQ(2)-g*sin(t1))/l1) + sum(4);
dQ(5) = dt2 + sum(5);
dQ(6) = (cos(t2)*dQ(2)-g*sin(t2))/l2 + sum(6);
end