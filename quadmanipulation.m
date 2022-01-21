clear all
close all
clc

syms theta1 theta2 theta1_dot theta2_dot a1 a2 m1 m2 g theta1_dd theta2_dd  

kf=9.8e-6;
km=1.6e-7;
m=2;
Ixx=0.0035;Iyy=0.0035;Izz=0.005;

l=0.225;
w1=100;w2=100;w3=100;w4=100;
u1=((w1^2)+(w2^2)+(w3^2)+(w4^2))*kf;

x_1 = 0; y_1 = 0;
x_2 = 0; y_2 = 0.15;

[Theta11, Theta12] = invk(x_1,y_1);
theta11 = real(Theta11), theta12 = real(Theta12)
[Theta21, Theta22] = invk(x_2,y_2);
theta21 = real(Theta21), theta22 = real(Theta22)

% tau1 =1; tau2 =1; 
theta1 = theta21; theta2 = theta22; 
theta1_dot = abs((theta21 - theta11)/2); theta2_dot = abs((theta22 - theta12)/2);
m1 = 0.1; m2 = 0.1; 
a1 = 0.10; a2 = 0.10; 
g = 9.81;

M = [m1*a1^2+m2*(a1^2+2*a1*a2*cos(theta2)+a2^2) m2*(a1*a2*cos(theta2)+a2^2);
     m2*(a1*a2*cos(theta2)+a2^2) 0]
S_arm = [-2*m2*a1*a2*sin(theta2)*theta2_dot -m2*a1*a2*sin(theta2)*theta2_dot;
     m2*a1*a2*sin(theta2)*theta1_dot 0]
G_arm = [(m1+m2)*a1*g*cos(theta1)+m2*g*a2*cos(theta1+theta2);
            m2*g*a2*cos(theta1+theta2)]
M_inv = inv(M);
MG = -M_inv*G_arm
A_arm = M_inv*S_arm;
k = -(2*a1*a2*theta2_dot*sin(theta2))/(a2^2 + a1*cos(theta2)*a2) - (a1*a2*theta1_dot*sin(theta2)*(a1^2*m1 + a1^2*m2 + a2^2*m2 + 2*a1*a2*m2*cos(theta2)))/(m2*(a2^2 + a1*cos(theta2)*a2)^2);
A=[0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
   0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
   0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 0 u1/m 0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 -u1/m 0 0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
   0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0 0 0 0 0 -(a1*theta1_dot*sin(theta2))/(a2 + a1*cos(theta2)) 0 0 0;
   0 0 0 0 0 0 0 0 0 0 0 0 -k (a1*theta2_dot*sin(theta2))/(a2 + a1*cos(theta2)) 0 0;
   0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
   0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0]

B=[0 0 0 0 0 0;
   0 0 0 0 0 0;
   0 0 0 0 0 0;
   0 0 0 0 0 0;
   0 0 0 0 0 0;
   1/m 0 0 0 0 0;
   0 0 0 0 0 0;
   0 0 0 0 0 0;
   0 0 0 0 0 0;
   0 01/Ixx 0 0 0 0;
   0 0 1/Iyy 0 0 0;
   0 0 0 1/Izz 0 0;
   0 0 0 0 0 1/(m2*(a2^2 + a1*cos(theta2)*a2));
   0 0 0 0 1/(m2*(a2^2 + a1*cos(theta2)*a2)) -(a1^2*m1 + a1^2*m2 + a2^2*m2 + 2*a1*a2*m2*cos(theta2))/(m2^2*(a2^2 + a1*cos(theta2)*a2)^2);
   0 0 0 0 0 0;
   0 0 0 0 0 0]
G = [0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    -(a2*g*cos(theta1 + theta2))/(a2^2 + a1*cos(theta2)*a2);
    (a2*g*cos(theta1 + theta2)*(a1^2*m1 + a1^2*m2 + a2^2*m2 + 2*a1*a2*m2*cos(theta2)))/(m2*(a2^2 + a1*cos(theta2)*a2)^2) - (a1*g*cos(theta1)*(m1 + m2) + a2*g*m2*cos(theta1 + theta2))/(m2*(a2^2 + a1*cos(theta2)*a2));
    0;
    0]

v= [500 200 1.71 500 200 3 10 10 10 0.25 10 1 100 100 1000 1000];
%x=[x y z x'y'z' phi theta zhi p    q    r]
Q=diag(v)
g=[1 0.001 0.001 0.001 1 1];
R=diag(g)
N=0;
[K,S,e] = lqr(A,B,Q,R,N);
J = [-a2*sin(theta1+theta2)-a1*sin(theta1) -a2*sin(theta1+theta2);
      a2*cos(theta1+theta2)+a1*cos(theta1) a2*cos(theta1+theta2)]
JT = transpose(J)
b_imp = [100 100]; 
B_imp = diag(b_imp)
k_imp = [100 1000];
K_imp = diag(k_imp)

function [theta1, theta2] = invk(x,y)
a1 = 0.10; 
a2 = 0.10;

d = sqrt(x^2 + y^2);
phi2 = acos((d^2 - a1^2 - a2^2)/(-2*a1*a2));
theta2 = pi - phi2;

phi1 = acos((a2^2 - a1^2 - d^2)/(-2*a1*d));

if(x==0)
    theta1 = pi/2 - phi1;
elseif(y<0 && x<0)
    theta1 = -pi + atan(y/x) - phi1;
else
    theta1 = atan(y/x) - phi1;
end

x1 = 0;
y1 = 0;

x2 = a1*cos(theta1);
y2 = a1*sin(theta1);

x3 = a2*cos(theta1+theta2) + x2;
y3 = a2*sin(theta1+theta2) + y2;

L1=sqrt((y2-y1)^2+(x2-x1)^2);
L2=sqrt((y3-y2)^2+(x3-x2)^2);
end