clear;
clc;

ts = 0.1;
tmax = 10;
N = 101;
T = 0:ts:(N-1)*ts;

s = tf('s');
G1 = 1/(s+1);
G2 = 10/(s+1);

Gd1 = c2d(G1, ts);
Gd2 = c2d(G2, ts);

B1 = conv(Gd1.num{1}(2:end), Gd2.den{1});
B2 = conv(Gd2.num{1}(2:end), Gd1.den{1});
A = conv(Gd1.den{1}, Gd2.den{1});
A = A(2:end);

nb = length(B1);
na = length(A);

Y11 = step(G1, T); 
Y12 = step(G2, T); 

Y = [Y11, Y12];

N = length(T);
Ny = 10;
N1 = 1;
Nu = 5;
lambda = 1;

[dmc] = dmc_setup(N, N1, Ny, Nu, lambda, Y)

k0 = 1+max([1+nb,1+na]);
kmax = tmax/ts + k0;

minLength = length([k0:kmax]) + k0 + Ny;

yr = 5*ones(minLength, 1);
yr(1:30) = 1;
yr(31:end) = 5;
yr(61:end) = 2;
u = zeros(minLength, 2);
y = zeros(minLength, 1);
dup = zeros(N-N1, 2);

for k = k0:kmax
    y(k) = B1*u(k-1:-1:k-nb, 1) + B2*u(k-1:-1:k-nb, 2) - A*y(k-1:-1:k-na);
    
    %du = dmc_controlaction_constrained(dmc, yr(k+[1:Ny]), dup, y(k), [u(k-1,1); u(k-1,2)], [-0.5 1; -0.5 1]);
    %du = dmc_controlaction_constrained(dmc, yr(k+[1:Ny]), dup, y(k), [], [], [-0.5 0.5; -0.5 0.5]);
    %du = dmc_controlaction_constrained(dmc, yr(k+[1:Ny]), dup, y(k), [u(k-1,1); u(k-1,2)], [-0.5 1; -0.5 1], [-0.5 0.5; -0.5 0.5]);
    %du = dmc_controlaction_constrained(dmc, yr(k+[1:Ny]), dup, y(k), [], [], [], [0 5.1]);
    du = dmc_controlaction(dmc, yr(k+[1:Ny]), dup, y(k));
    
    u(k, 1) = u(k-1, 1) + du(1);
    u(k, 2) = u(k-1, 2) + du(2);
    
    dup(2:end, :) = dup(1:end-1, :);
    dup(1, :) = du';
end

figure;
subplot(2,1,1);
stairs(([k0:kmax]' - k0)*ts, yr(k0:kmax), 'k--');
hold all;
stairs(([k0:kmax]' - k0)*ts, y(k0:kmax), 'k');
subplot(2,1,2);
stairs(([k0:kmax]' - k0)*ts, u(k0:kmax, :));