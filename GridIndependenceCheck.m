clear; clc

% Grid resolutions
Nx = [100 140 196 276 388];
Ny = [16 24 36 52 72];
N = [1717 3525 7289 14681 28397];

% Parameters for design 3
Nx1 = Nx/2;
Ny1 = Ny/2;
k1 = 100; k3 = 100; k4 = 100;
k2 = 1000;

% Computing q_f for every grid
q_f = zeros(1, 5);
for n = 1:5
    [T, q_f(1, n)] = NumericalSolution(Nx(n), Nx1(n), Ny(n), Ny1(n), k1, k2, k3, k4);
end

% Computing error relative to the 5th grid 
Er_N = abs((q_f(1, 1:4) - q_f(1, 5))/q_f(1, 5)) * 100;

% Plotting the result
loglog(N(1:4), Er_N)
grid on
title("log-log Plot of Er vs. N")
xlabel("N(total node number)")
ylabel("Er(%)")