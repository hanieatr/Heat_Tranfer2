function [T, q_f] = AnalyticalSolution(Nx)
t = 10e-3; % Thickness
L = 100e-3; % Length
w = 1; % Width
h = 200; % Convection coefficient
T_inf = 300; % Environment temperature
Tb = 350; % Base temperature
k = 100; % Conductivity
x = linspace(0, L, Nx + 1); % Analytical solution will be compared to a numerical solution with Nx+1 node(Nx=400)

P = 2*w + 2*t; % Fin perimeter
A_c = w*t; % Fin cross-sectional area

% Other parameters for analytical solution
m = sqrt(h*P/(k*A_c));
theta_b = Tb - T_inf;
M = sqrt(h*P*k*A_c) * theta_b;

% Temperature dist
T = T_inf + theta_b * (cosh(m*(L-x)) + (h/(m*k))*sinh(m*(L-x)))./(cosh(m*L) + (h/(m*k))*sinh(m*L));
% Heat transfer rate
q_f = M * (sinh(m*L) + (h/(m*k))*cosh(m*L))/(cosh(m*L) + (h/(m*k))*sinh(m*L));
end