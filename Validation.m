clear; clc

% Parameters for validation
L = 100e-3; % Length
Nx = 400; % Total number of sections in x-direction
Nx1 = 200; % Number of sections in x-direction(First part)
Ny = 20; % Total number of sections in y-direction
Ny1 = 10; % Number of sections in y-direction(First part)
k = 100; % Conductivity
x = linspace(0, L, Nx+1);

% Numerical solution
[T_numerical, q_numerical] = NumericalSolution(Nx, Nx1, Ny, Ny1, k, k, k, k);

% Analytical solution
[T_analytical, q_analytical] = AnalyticalSolution(Nx);
T_y1 = T_numerical(:, 1)'; % y = 0 (y is measured from the centerline of the fin)
T_y2 = T_numerical(:, 11)'; % y = t/4 (t is the thickness of the whole fin)
T_y3 = T_numerical(:, end)'; % y = t/2

% Plotting the result
plot(x, T_analytical, x, T_y1, x, T_y2, x, T_y3)
xlabel("x(m)")
ylabel("T(K)")
legend("Analytical", "Numerical(y = 0)", "Numerical(y = t/4)", "Numerical(y = t/2)")