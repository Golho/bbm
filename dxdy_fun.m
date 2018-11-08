function [dxdy] = dxdy_fun(y, X)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
A = 0.0002359;
B = -0.6355;
C = -0.6389;
D = 0.6953;

h = 0.3;
L = 6;
rho_0 = 1;
mu = -1;
mu_2 = -1;
theta = deg2rad(30);
g = 9.81;
g_x = cos(theta)*g;

rho = @(y) rho_0*(A*exp(y*L/h) + B*exp(-y*L/h) + C*y/h + D);
drhody = @(y) rho_0*(L/h*A*exp(y*L/h) - L/h*B*exp(-y*L/h) + C/h);
temp2 = @(y) 2*rho_0^2*L^2/h^3*(L*A^2*exp(2*L/h*y) + A*C*exp(L/h*y) - L*B^2*exp(-2*L/h*y) + B*C*exp(-L/h*y));

f = @(y) mu + mu_2/2*drhody(y)^2;
g = @(y) mu_2/2 * temp2(y);
h = @(y) -rho(y)*g_x;
dxdy(1) = X(2);
dxdy(2) = (h(y) - g(y)*X(2))/f(y);
end

