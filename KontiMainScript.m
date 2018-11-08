% Assignment in Continuum Mechanics 2018
clear all;
clc;
close all;
%% Exercise 1)

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
i = 1;
for L=4:0.5:8
    syms A B C D
    
    % Symbolically express all the boundary conditions equations
    ekv1 = A + B + D == v_0y;
    ekv2 = A*exp(0.3*L) + B*exp(-0.3*L) + 0.3*C + D == v_03y;
    ekv3 = A*exp(L) + B*exp(-L) + C + D == v_1y;
    ekv4 = L*(A*exp(0.3*L)-B*exp(-0.3*L)) + C == 0;

    % Symbolically solve the boundary conditions equations for each L
    s = solve([ekv1,ekv2,ekv3,ekv4],[A B C D]);

    % Plot the solution
    f = @(x) (s.A*exp(x*L)+s.B*exp(-x*L)+s.C*x+s.D);
    plotvector(i) = fplot(f,[0 1]);
    i = i + 1;
end
axis([0 1 0 0.8])
legend(plotvector, '4', '4.5', '5', '5.5', '6', '6.5', '7', '7.5', '8')
hold off

figure(2)
hold on
title('Solids fraction \nu vs nondimensional depth ')
syms A B C D
xlabel('y/h');
ylabel('\nu');

% Choose L = 6 and solve the equation system
L = 6;
ekv1 = A + B + D == v_0y;
ekv2 = A*exp(0.3*L) + B*exp(-0.3*L) + 0.3*C + D == v_03y;
ekv3 = A*exp(L) + B*exp(-L) + C + D == v_1y;
ekv4 = L*(A*exp(0.3*L)-B*exp(-0.3*L)) + C == 0;

s2 = solve([ekv1,ekv2,ekv3,ekv4],[A B C D]);
nu = @(y) (s2.A*exp(y*L)+s2.B*exp(-y*L)+s2.C*y+s2.D);
fplot(nu,[0 1])
hold off

%% Exercise 2

% Integration constant
k = h/L*(s2.B-s2.A);
rho_0 = 1;
theta = 30;
g = 9.81;
b = g*[sind(theta), cosd(theta), 0]; % [b_x, b_y, b_z]
my = -1;
g_y = cosd(theta)*g;
g_x = sind(theta)*g;
p_1 = @(y) (rho_0*b(2)*(s2.A*h/L*exp(y/h*L)-s2.B*h/L*exp(-y/h*L)+s2.C*y.^2/(2*h)+s2.D*y + k));
figure(3)
hold on
title('Gauge pressure vs flow depth')
xlabel('y [m]');
ylabel('p [Pa (Gauge)]');
fplot(p_1,[0 h])
hold off
% Integration constant
l = -(h/L)^2*(s2.A*exp(L)+s2.B*exp(-L)+s2.C*L^2/6+s2.D*L^2/2)-k*h;
xdot_1 = @(x) -b(1)*rho_0/my * (s2.A*(h/L)^2*exp(x/h*L)+s2.B*(h/L)^2*exp(-x/h*L)+s2.C*x.^3/(6*h)+s2.D*x.^2/2+k*x+l);
figure(4)
hold on
title('Flow velocity vs flow depth')
xlabel('y [m]');
ylabel('velocity [m/s]');
fplot(xdot_1,[0 h])
hold off
%% Exercise 3
A = s2.A;
B = s2.B;
C = s2.C;
D = s2.D;
mu_1 = 0.005;
p_2 = @(y) rho_0^2*(L/h)^2*(A^2*(exp(2*L*y/h)-1) + B^2*(exp(-2*L*y/h)-1) + 2/L*A*C*(exp(L*y/h)-1) - 2/L*B*C*(exp(-L*y/h)-1));
p_const2 = @(y) p_1(y) + mu_1*p_2(y);
figure
fplot(p_const2,[0 h])
title('Gauge pressure vs flow depth')
xlabel('y [m]');
ylabel('p [Pa (Gauge)]');
%% Solve differential equation
clear h
h = 0.3;
options = bvpset('RelTol', 1e-5);
solinit = bvpinit([0,h],[0,0.12]);
 
sol = bvp4c(@dxdy_fun, @bcs, solinit, options);
figure
plot(sol.x,sol.y(1,:));
title('Flow velocity vs flow depth')
xlabel('y [m]');
ylabel('velocity [m/s]');

function [dxdy] = dxdy_fun(y, X)
%dxdu_fun Return the derivative of the flow velocity for constitutive
%relation 2
A = 0.0002359;
B = -0.6355;
C = -0.6389;
D = 0.6953;

h = 0.3;
L = 6;
rho_0 = 1;
mu = -1;
mu_2 = -1;
theta = 30;
g = 9.81;
g_x = sind(theta)*g;

rho = @(y) rho_0*(A*exp(y*L/h) + B*exp(-y*L/h) + C*y/h + D);
drhody = @(y) rho_0*(L/h*A*exp(y*L/h) - L/h*B*exp(-y*L/h) + C/h);
temp2 = @(y) 2*rho_0^2*L^2/h^3*(L*A^2*exp(2*L/h*y) + A*C*exp(L/h*y) - L*B^2*exp(-2*L/h*y) + B*C*exp(-L/h*y));

f = @(y) mu + mu_2/2*drhody(y)^2;
g = @(y) mu_2/2 * temp2(y);
h = @(y) -rho(y)*g_x;
dxdy(1) = X(2);
dxdy(2) = (h(y) - g(y)*X(2))/f(y);
end

function [res] = bcs(xa,xb)
%bcs Return the boundary conditions for the flow for constitutive relation
% 2
res = [xa(2);
       xb(1)];
end
