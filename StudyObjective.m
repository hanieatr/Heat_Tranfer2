clear; clc

t = 5e-3; % Thickness
L = 100e-3; % Length
Nx = 100; % Total number of sections in x-direction
Ny = 16; % Total number of sections in y-direction
k_a = 100; % Conductivity of Aluminum alloy
k_CNF = 1000; % Conductivity of Carbon NanoFiber

% Parameters for different designs
Nx1 = [Nx/2, Nx/2, Nx/2, Nx/2, Nx/2, Nx/2, Nx/4, 3*Nx/4];
Ny1 = [3*Ny/4, Ny/4, Ny/2, Ny/2, Ny/2, Ny/2, Ny/2, Ny/2];
k1 = [k_a, k_CNF, k_a, k_CNF, k_a, k_a, k_CNF, k_a];
k2 = [k_CNF, k_a, k_CNF, k_a, k_a, k_a, k_CNF, k_a];
k3 = [k_a, k_CNF, k_a, k_a, k_a, k_CNF, k_a, k_CNF];
k4 = [k_CNF, k_a, k_a, k_a, k_CNF, k_a, k_a, k_CNF];

% Temperature and heat transfer rate for base design
[T_basedesign, q_f0] = NumericalSolution(Nx, Nx/2, Ny, Ny/2, 100, 100, 100, 100);

% Temperature and heat transfer rate for all designs
T = zeros(Nx+1, Ny+1, 8);
q_f = zeros(1, 8);
for n=1:8
    [T(:, :, n), q_f(1, n)] = NumericalSolution(Nx, Nx1(n), Ny, Ny1(n), k1(n), k2(n), k3(n), k4(n));
end

% Computing the ratio of the heat transfer rates
ratios = q_f/q_f0;
[max_ratio, argmax_ratio] = max(ratios);

% Plotting the result
x = categorical(["Design1", "Design2", "Design3", "Design4", "Design5", "Design6", "Design7", "Design8"]);
b = bar(x, ratios, "FaceColor", "flat");
b.CData(argmax_ratio,:) = [0 0.8 0.8];
ylabel("q_f/q_f_0")
title("q_f/q_f_0 Ratio for Different Designs")