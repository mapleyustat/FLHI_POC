clear;
clc;

% Simulation time parameters
ts = 0.01;
tmax = 5;

% Step model parameters
N = 501;
T = 0:ts:(N-1)*ts;

% DMC parameters
N = length(T);
Ny = N;
N1 = 1;
Nu = 10;
lambda = 10;

% MIMO step model
s = tf('s');
G1 = 1/(s+1);
G2 = 10/(s+1);
G3 = 1/(s+1);
G4 = 10/(s+10);

Gd1 = c2d(G1, ts);
Gd2 = c2d(G2, ts);
Gd3 = c2d(G3, ts);
Gd4 = c2d(G4, ts);

B1 = conv(Gd1.num{1}(2:end), Gd2.den{1});
B2 = conv(Gd2.num{1}(2:end), Gd1.den{1});
A1 = conv(Gd1.den{1}, Gd2.den{1});
A1 = A1(2:end);

B3 = conv(Gd3.num{1}(2:end), Gd4.den{1});
B4 = conv(Gd4.num{1}(2:end), Gd3.den{1});
A2 = conv(Gd3.den{1}, Gd4.den{1});
A2 = A2(2:end);

nb = length(B1);
na = length(A1);

Y11 = step(G1, T); 
Y12 = step(G2, T); 
Y21 = step(G3, T); 
Y22 = step(G4, T); 

input_count = 2;
output_count = 2;
Y = zeros(length(T), input_count, output_count);
% Y(:,:,1) = [Y11, Y12];
% Y(:,:,2) = [Y21, Y22];
Y(:,1,1) = Y11;
Y(:,2,1) = Y12;
Y(:,1,2) = Y21;
Y(:,2,2) = Y22;

% Create a DMC controller for the step models
[dmc] = dmc_setup(N, N1, Ny, Nu, lambda, Y)

% Setup additional simulation variables
k0 = 1+max([1+nb,1+na]);
kmax = tmax/ts + k0;

minLength = length([k0:kmax]) + k0 + Ny;

yr = 5*ones(minLength, output_count);
% yr(1:300, :) = 1;
% yr(301:end, :) = 5;
% yr(601:end, :) = 2;
u = zeros(minLength, input_count);
y = zeros(minLength, output_count);
dup = zeros(N-N1, input_count);

f = [0 0];

% Main simulation loop
for k = k0:kmax
    y(k, 1) = B1*u(k-1:-1:k-nb, 1) + B2*u(k-1:-1:k-nb, 2) - A1*y(k-1:-1:k-na, 1);
    y(k, 2) = B3*u(k-1:-1:k-nb, 1) + B4*u(k-1:-1:k-nb, 2) - A2*y(k-1:-1:k-na, 2);
    
    %du = dmc_controlaction_constrained(dmc, yr(k+[1:Ny-N1+1], :), dup, y(k, :), [u(k-1,1); u(k-1,2)], [-0.5 5; -0.5 5]);
    %du = dmc_controlaction_constrained(dmc, yr(k+[1:Ny-N1+1], :), dup, y(k, :), [], [], [-0.5 0.5; -0.5 0.5]);
    %du = dmc_controlaction_constrained(dmc, yr(k+[1:Ny-N1+1], :), dup, y(k, :), [u(k-1,1); u(k-1,2)], [-0.5 5; -0.5 5], [-0.5 0.5; -0.5 0.5]);
    %du = dmc_controlaction_constrained(dmc, yr(k+[1:Ny-N1+1], :), dup, y(k,:), [], [], [], [-10 6; -10 6]);
    du = dmc_controlaction(dmc, yr(k+[1:Ny-N1+1], :), dup, y(k,:));
    
    u(k, 1) = u(k-1, 1) + du(1);
    u(k, 2) = u(k-1, 2) + du(2);
    
    dup(2:end, :) = dup(1:end-1, :);
    dup(1, :) = du';
end

figure;
subplot(2,1,1);
stairs(([k0:kmax]' - k0)*ts, yr(k0:kmax,:), 'k--');
hold all;
ylabel('Saída (mag)');
xlabel('Tempo (t)');
stairs(([k0:kmax]' - k0)*ts, y(k0:kmax,:), 'LineWidth', 2);
subplot(2,1,2);
stairs(([k0:kmax]' - k0)*ts, u(k0:kmax, :), 'LineWidth', 2);
ylabel('Controle (mag)');
xlabel('Tempo (t)');