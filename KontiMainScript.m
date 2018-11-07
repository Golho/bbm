       
              % Assignment in Continuum Mechanics 2018
clear all;
clc;
close all;
%% excercise 1)

% From curve 1 (solid line)
v_0y = 0.06;
v_03y = 0.4;
v_1y = 0.15;
h = 0.3;

figure(1)
hold on
xlabel('v');
ylabel('y/h');
plotvector = zeros(1,9);
%L = 1;
i = 1;
for L=4:0.5:8
    syms A B C D

        ekv1 = A + B + D == v_0y;
        ekv2 = A*exp(0.3*L) + B*exp(-0.3*L) + 0.3*C + D == v_03y;
        ekv3 = A*exp(L) + B*exp(-L) + C + D == v_1y;
        ekv4 = L*(A*exp(0.3*L)-B*exp(-0.3*L)) + C == 0;

        s = solve([ekv1,ekv2,ekv3,ekv4],[A B C D]);
        f = @(x) (s.A*exp(x*L)+s.B*exp(-x*L)+s.C*x+s.D);
        plotvector(i) = fplot(f,[0 1]);
        i = i + 1;
    
end
axis([0 1 0 0.8])
legend(plotvector, '4', '4.5', '5', '5.5', '6', '6.5', '7', '7.5', '8')
hold off
%legend ('1','2',... - make for loop to assemble string
% L = 6
figure(2)
hold on
title('Solids fraction \nu vs nondimensional depth ')
syms A B C D
xlabel('y/h');
ylabel('\nu');

L = 6;
ekv1 = A + B + D == v_0y;
ekv2 = A*exp(0.3*L) + B*exp(-0.3*L) + 0.3*C + D == v_03y;
ekv3 = A*exp(L) + B*exp(-L) + C + D == v_1y;
ekv4 = L*(A*exp(0.3*L)-B*exp(-0.3*L)) + C == 0;

s2 = solve([ekv1,ekv2,ekv3,ekv4],[A B C D]);
f_2 = @(x) (s2.A*exp(x*L)+s2.B*exp(-x*L)+s2.C*x+s2.D);
fplot(f_2,[0 1])
hold off

%% Excercise 2

k = h/L*(s2.B-s2.A);
rho_0 = 1;
theta = 30;
g = 9.81;
my = 1;
g_y = cosd(theta)*g;
g_x = sind(theta)*g;
f_p = @(x) (rho_0*g_x*(s2.A*h/L*exp(x/h*L)-s2.B*h/L*exp(-x/h*L)+s2.C*x.^2/(2*h)+s2.D*x+k));
figure(3)
hold on
title('Gauge pressure vs flow depth')
xlabel('y [m]');
ylabel('p [Pa (Gauge)]');
fplot(f_p,[0 h])
hold off
l = -(h/L)^2*(s2.A*exp(L)+s2.B*exp(-L)+s2.C*L^2/6+s2.D*L^2/2)-k*h;
f_xdot = @(x) (-g_y/my*(s2.A*(h/L)^2*exp(x/h*L)+s2.B*(h/L)^2*exp(-x/h*L)+s2.C*x.^3/(6*h)+s2.D*x.^2/2+k*x+l));
figure(4)
hold on
title('Flow velocity vs flow depth')
xlabel('y [m]');
ylabel('velocity [m/s]');
fplot(f_xdot,[0 h])
hold off
