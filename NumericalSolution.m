function [Tkp1, q_f] = NumericalSolution(Nx, Nx1, Ny, Ny1, k1, k2, k3, k4)
% Other inputs
t = 5e-3; % Thickness
L = 100e-3; % Length
h = 200; % Convection coefficient
T_inf = 300; % Environment temperature
Tb = 350; % Base temperature

% Numerical Parameters Setup
dx = L/Nx; 
dy = t/Ny;
R = dx^2/dy^2;

% First guess as T_b(by defining the first guess like this, we won't need
% to define the isothermal boundary condition again
Tk = Tb*ones(Nx+1, Ny+1); % Temperature distribution at iteration k
Tkp1 = Tb*ones(Nx+1, Ny+1); % Temperature distribution at iteration k+1

% Heat transfer rate
q_f = 0;

% Main loop(solving with Gauss-Siedel method)
while true
    % interior
    for n = 2:Ny
        for m = 2:Nx
            if n~=Ny1+1 && m~=Nx1+1
                Tkp1(m,n) = 1/(2*(1+R)) * (Tkp1(m-1,n) + Tkp1(m+1,n) + ...
                    R * (Tkp1(m,n+1) + Tkp1(m,n-1)));
            end
        end
    end

    % eq1
    n = Ny + 1;
    for m = 2:Nx
        if m <= Nx1
            Tkp1(m,n) = 1/(dx^2/dy*h + k2 + k2*R) *...
                (dx^2/dy*h*T_inf + 0.5*k2*(Tkp1(m+1,n) + Tkp1(m-1,n)) + k2*R*Tkp1(m,n-1));
        elseif m > Nx1 + 1
            Tkp1(m,n) = 1/(dx^2/dy*h + k4 + k4*R) *...
                (dx^2/dy*h*T_inf + 0.5*k4*(Tkp1(m+1,n) + Tkp1(m-1,n)) + k4*R*Tkp1(m,n-1));
        end
    end

    % eq2
    n = Ny + 1;
    m = Nx1 + 1;
    Tkp1(m,n) = 1/(dx^2/dy*h + 0.5*(k2 + k4)*(1 + R)) *...
        (dx^2/dy*h*T_inf + 0.5*(k4*Tkp1(m+1,n) + k2*Tkp1(m-1,n)) + 0.5*R*(k2 + k4)*Tkp1(m,n-1));

    % eq3
    n = Ny + 1;
    m = Nx + 1;
    Tkp1(m,n) = 1/((dx^2/dy + dx)*h + k4*(1 + R)) *...
        ((dx^2/dy + dx)*h*T_inf + k4*R*Tkp1(m,n-1) + k4*Tkp1(m-1,n));

    % eq4
    n = Ny1 + 1;
    m = Nx + 1;
    Tkp1(m,n) = 1/(dx*h + 0.5*(k3 + k4)*(1 + R)) *...
        (dx*h*T_inf + 0.5*(k3+k4)*Tkp1(m-1,n) + 0.5*R*(k4*Tkp1(m,n+1) + k3*Tkp1(m,n-1)));

    % eq5
    n = Ny1 + 1;
    m = Nx1 + 1;
    Tkp1(m,n) = 1/((k1+k2+k3+k4)*(1 + R)) *...
        ((k3+k4)*Tkp1(m+1,n) + (k1+k2)*Tkp1(m-1,n) + R*(k2+k4)*Tkp1(m,n+1) + R*(k1+k3)*Tkp1(m,n-1));

    % eq6
    n = Ny1 + 1;
    for m = 2:Nx
        if m <= Nx1
            Tkp1(m,n) = 1/((k1+k2)*(1 + R)) *...
                (0.5*(k1+k2)*(Tkp1(m+1,n) + Tkp1(m-1,n)) + k2*R*Tkp1(m,n+1) + k1*R*Tkp1(m,n-1));
        elseif m > Nx1 + 1
            Tkp1(m,n) = 1/((k3+k4)*(1 + R)) *...
                (0.5*(k3+k4)*(Tkp1(m+1,n) + Tkp1(m-1,n)) + k4*R*Tkp1(m,n+1) + k3*R*Tkp1(m,n-1));
        end
    end

    % eq7
    m = Nx + 1;
    for n = 2:Ny
        if n <= Ny1
            Tkp1(m,n) = 1/(dx*h + (1 + R)*k3) *...
                (dx*h*T_inf + k3*Tkp1(m-1,n) + 0.5*R*k3*(Tkp1(m,n+1) + Tkp1(m,n-1)));
        elseif n > Ny1 + 1
            Tkp1(m,n) = 1/(dx*h + (1 + R)*k4) *...
                (dx*h*T_inf + k4*Tkp1(m-1,n) + 0.5*R*k4*(Tkp1(m,n+1) + Tkp1(m,n-1)));
        end
    end

    % eq8
    n = 1;
    m = Nx + 1;
    Tkp1(m,n) = 1/(dx*h + (1 + R)*k3) *...
                (dx*h*T_inf + k3*Tkp1(m-1,n) + R*k3*Tkp1(m,n+1));

    % eq9
    m = Nx1 + 1;
    for n = 2:Ny
        if n <= Ny1
            Tkp1(m,n) = 1/((k1+k3)*(1 + R)) *...
                (k3*Tkp1(m+1,n) + k1*Tkp1(m-1,n) + 0.5*R*(k1+k3)*(Tkp1(m,n+1) + Tkp1(m,n-1)));
        elseif n > Ny1 + 1
            Tkp1(m,n) = 1/((k2+k4)*(1 + R)) *...
                (k4*Tkp1(m+1,n) + k2*Tkp1(m-1,n) + 0.5*R*(k2+k4)*(Tkp1(m,n+1) + Tkp1(m,n-1)));
        end
    end

    % eq10
    n = 1;
    m = Nx1 + 1;
    Tkp1(m,n) = 1/((k1+k3)*(1 + R)) *...
                (k3*Tkp1(m+1,n) + k1*Tkp1(m-1,n) + R*(k1+k3)*Tkp1(m,n+1));

    % eq11
    n = 1;
    for m = 2:Nx
        if m <= Nx1
            Tkp1(m,n) = 1/(k1*(1 + R)) *...
                (0.5*k1*(Tkp1(m+1,n) + Tkp1(m-1,n)) + k1*R*Tkp1(m,n+1));
        elseif m > Nx1 + 1
            Tkp1(m,n) = 1/(k3*(1 + R)) *...
                (0.5*k3*(Tkp1(m+1,n) + Tkp1(m-1,n)) + k3*R*Tkp1(m,n+1));
        end
    end

    
    % Computing the 2nd norm of error
    er = (Tkp1-Tk).^2;
    er2 = sqrt(sum(er, 'all'));

    % Checking for convergence
    if er2 <= 1e-5
        break
    end	
	
    % Preparing for the next iteration
	Tk = Tkp1; 
end


% Computing q_f
for n = 1:Ny+1
    if n == 1
        q_f = q_f + k1 * dy/2 * (Tkp1(1,n)-Tkp1(2,n))/dx;
    elseif n <= Ny1
        q_f = q_f + k1 * dy * (Tkp1(1,n)-Tkp1(2,n))/dx;
    elseif n == Ny1 + 1
        q_f = q_f + (k1 + k2) * dy/2 * (Tkp1(1,n)-Tkp1(2,n))/dx;
    elseif n > Ny1 + 1 && n < Ny + 1
        q_f = q_f + k2 * dy * (Tkp1(1,n)-Tkp1(2,n))/dx;
    else
        q_f = q_f + k2 * dy/2 * (Tkp1(1,n)-Tkp1(2,n))/dx;
    end
end
q_f = q_f + h * dx/2 * (Tkp1(1,Ny+1) - T_inf); % Convection term
q_f = q_f * 2; % Computing for the whole fin
end