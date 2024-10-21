clear; clc

t = 5e-3; % Thickness
L = 100e-3; % Length
Nx = 100; % Total number of sections in x-direction
Ny = 16; % Total number of sections in y-direction

% Parameters for design 3
Nx1 = Nx/2;
Ny1 = Ny/2;
k1 = 100; k3 = 100; k4 = 100;
k2 = 1000;

% Temperature and heat transfer rate for base design
[T_basedesign, q_f0] = NumericalSolution(Nx, Nx1, Ny, Ny1, 100, 100, 100, 100);

% Temperature and heat transfer rate for design 3
[T_numerical, q_f] = NumericalSolution(Nx, Nx1, Ny, Ny1, k1, k2, k3, k4);

% Computing the ratio of the heat transfer rates
ratio = q_f/q_f0;


% Plotting the result
x=linspace(0, L, Nx+1);
y=linspace(0, t, Ny+1);
[X,Y] = meshgrid(x, y);

% Design 3 contour
figure(1)
[c, h] = contourf(X, Y, T_numerical');
clabel(c, h);
hbar=colorbar;
ylabel(hbar, 'T(K)');
colormap(jet);
xlabel('x(m)');
ylabel('y(m)');
title('Temperature Distribution(Design 3)');

% Base design contour
figure(2)
[c, h] = contourf(X, Y, T_basedesign');
clabel(c, h);
hbar=colorbar;
ylabel(hbar, 'T(K)');
colormap(jet);
xlabel('x(m)');
ylabel('y(m)');
title('Temperature Distribution(Base Design)');