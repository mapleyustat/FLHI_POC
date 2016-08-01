clear;
clc;

s = tf('s');
G = 1/(s+1);

ts = 0.1;
tmax = 10;
N = 101;
T = 0:ts:(N-1)*ts;

Gd = c2d(G, ts);

B = Gd.num{1}(2:end);
A = Gd.den{1}(2:end);

nb = length(B);
na = length(A);

Y11 = step(G, T); 

N = length(T);
Ny = 5;
N1 = 1;
Nu = 1;
lambda = 1;

[dmc] = dmc_setup(N, N1, Ny, Nu, lambda, Y11)

k0 = 1+max([1+nb,1+na]);
kmax = tmax/ts + k0;

minLength = length([k0:kmax]) + k0 + Ny;

yr = 5*ones(minLength, 1);
yr(1:30) = 1;
yr(31:end) = 5;
yr(61:end) = 2;
u = zeros(minLength, 1);
y = zeros(minLength, 1);
dup = zeros(N-N1, 1);

for k = k0:kmax
    y(k) = B*u(k-1:-1:k-nb) - A*y(k-1:-1:k-na);
    
    %du = dmc_controlaction_constrained(dmc, yr(k+[1:Ny]), dup, y(k), u(k-1), [0 7]);
    %du = dmc_controlaction_constrained(dmc, yr(k+[1:Ny]), dup, y(k), [], [], [-1 1]);
    %du = dmc_controlaction_constrained(dmc, yr(k+[1:Ny]), dup, y(k), u(k-1), [0 7], [-1 1]);
    %du = dmc_controlaction_constrained(dmc, yr(k+[1:Ny]), dup, y(k), [], [], [], [0 5]);
    du = dmc_controlaction_constrained(dmc, yr(k+[1:Ny]), dup, y(k), u(k-1), [0 7], [], [0 5]);
    %du = dmc_controlaction(dmc, yr(k+[1:Ny]), dup, y(k));
    
    u(k) = u(k-1) + du;
    
    dup(2:end) = dup(1:end-1);
    dup(1) = du;
end

figure;
subplot(2,1,1);
stairs(([k0:kmax]' - k0)*ts, yr(k0:kmax), 'k--');
hold all;
stairs(([k0:kmax]' - k0)*ts, y(k0:kmax), 'k');
subplot(2,1,2);
stairs(([k0:kmax]' - k0)*ts, u(k0:kmax));