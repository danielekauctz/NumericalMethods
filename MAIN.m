%% EXAMPLES
clear all; close all;
clc
%% EXAMPLE 1
% Problem data
m = 10;      
k = 40;      
z = 0.01;
x0 = 0.5;    
v0 = 0;

% Calculus
wn = sqrt(k/m);
c = 2*m*z*wn;

DOF = 1;
n = 6;
f = zeros(DOF,n);
dt = 0.2;
t0 = 0;
tf = n*dt;
t = linspace(t0,tf,n);


[x1,v1,a1] = const_acc_method(m,k,c,x0,v0,f,DOF,dt,n);
[x2,v2,a2] = finite_diff_method(m,k,c,x0,v0,f,DOF,dt,n);
[x3,v3,a3] = newmark_method(m,k,c,x0,v0,f,DOF,dt,n);

figure(1)
plot(t,x1)
xlabel('\fontsize{11}\bf Time (s)');
ylabel('\fontsize{11}\bf Displacement (m)');
title('\fontsize{13}\bf Constant Acceleration Method')

figure(2)
plot(t,x2)
xlabel('\fontsize{11}\bf Time (s)');
ylabel('\fontsize{11}\bf Displacement (m)');
title('\fontsize{13}\bf Central Difference Method')

figure(3)
plot(t,x3)
xlabel('\fontsize{11}\bf Time (s)');
ylabel('\fontsize{11}\bf Displacement (m)');
title('\fontsize{13}\bf Newmark''s Method')

%% EXAMPLE 2
% Problem data
m1 = 8;
m2 = 4;
M = [m1   0;
      0  m2];
  
k1 = 60;     
k2 = 30;
K = [k1+k2 -k2;
       -k2  k2]; 
   
C = [0 0;
     0 0];
 
x10 = 0.5;  
x20 = 0;
x0 = [x10 x20]';

v10 = 0;
v20 = 0;
v0 = [v10 v20]';

% Calculus
DOF = 2;
n = 6;
F = zeros(DOF,n);
dt = 0.2;
t0 = 0;
tf = n*dt;
t = linspace(t0,tf,n);

[x4,v4,a4] = const_acc_method(M,K,C,x0,v0,F,DOF,dt,n);
[x5,v5,a5] = finite_diff_method(M,K,C,x0,v0,F,DOF,dt,n);
[x6,v6,a6] = newmark_method(M,K,C,x0,v0,F,DOF,dt,n);

figure(4)
plot(t,x4(1,:),'b'); hold on;
plot(t,x4(2,:),'r'); hold off;
xlabel('\fontsize{11}\bf Time (s)');
ylabel('\fontsize{11}\bf Displacement (m)');
legend('Mass 1','Mass 2');
title('\fontsize{13}\bf Constant Acceleration Method')

figure(5)
plot(t,x5(1,:),'b'); hold on;
plot(t,x5(2,:),'r'); hold off;
xlabel('\fontsize{11}\bf Time (s)');
ylabel('\fontsize{11}\bf Displacement (m)');
legend('Mass 1','Mass 2');
title('\fontsize{13}\bf Central Difference Method')

figure(6)
plot(t,x6(1,:),'b'); hold on;
plot(t,x6(2,:),'r'); hold off;
xlabel('\fontsize{11}\bf Time (s)');
ylabel('\fontsize{11}\bf Displacement (m)');
legend('Mass 1','Mass 2');
title('\fontsize{13}\bf Newmark''s Method')
