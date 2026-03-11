%% Heat Conduction Lab Part 2 - Task 1
% Author: Sam Kleiner
clc;
clear;
close all;

% Givens (from pt 1)
x0_toTC1_in = 1 + (3/8); % [in]
x0_toTC1 = x0_toTC1_in * 0.0254; % [m]

dx_in = 1.0; % [in]
dx = dx_in * 0.0254; % [m]
nTC = 8;
x = x0_toTC1 + (0:nTC-1)*dx;

L = x(end);
x_eval = x(end);

% Aluminum properties:
k = 130;
rho = 2810;
cp = 960;

alpha = k/(rho*cp);

T0 = 17.06;
H = 91.08;

times = [1 1000];
Nmax = 10;
modes = 0:Nmax;
T_modes = zeros(Nmax+1,length(times));

%% Convergence
for j = 1:length(times)
    t = times(j);
    for N = 0:Nmax
        T = T0 + H*x_eval; % steady-state
        for n = 1:N  % transient
            lambda = (2*n-1)*pi/(2*L);
            bn = (-8*H*L)/((2*n-1)^2*pi^2) * (-1)^(n-1);
            T = T + bn*sin(lambda*x_eval)*exp(-alpha*lambda^2*t);
        end
        T_modes(N+1,j) = T;
    end
end

%% Plot
figure();
plot(modes,T_modes(:,1),'o-','LineWidth',2);
hold on;
plot(modes,T_modes(:,2),'s-','LineWidth',2);
grid on;
xlabel('Number of modes');
ylabel('Temperature at thermocouple (C)');
legend('t = 1 s','t = 1000 s','Location','best');
title('Analytical Solution for Convergence');

%% Fourier Number
Fo1 = alpha*times(1)/L^2;
Fo2 = alpha*times(2)/L^2;

fprintf('t = 1 s: %.4f\n',Fo1)
fprintf('t = 1000 s: %.4f\n',Fo2)