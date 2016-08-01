clearvars -except FIGURE_NAME INTERPOLANT_TYPE OUTPUT_YMAX OUTPUT_YMIN CONTROL_YMAX CONTROL_YMIN LIMIT_AXES IGNORE_NONLINEARITY

disp('====================');
fprintf('=  DMC SISO %s\n', FIGURE_NAME);
disp('====================');

ts = 1;
tmax = 100;
N = 26;
T = 0:ts:(N-1)*ts;

%Gd = c2d(G, ts);
Gd = filt([0, 0.243], [1, -0.757], ts);

uss = [-0.15:0.01:0.15];
vss = 1.04.*uss - 14.11.*uss.^2 - 16.72.*uss.^3 + 562.75.*uss.^4;

pointVector = [];

for i = 1:length(uss)
    point.xi = uss(i);
    point.xo = vss(i);
    
    pointVector = [pointVector; point];
end

[interpolant] = flhiInterpolant(pointVector, INTERPOLANT_TYPE);

B = Gd.num{1}(2:end);
A = Gd.den{1}(2:end);

nb = length(B) - 1;
na = length(A);

Y11 = step(Gd, T);
Y11 = Y11(2:end); % skip unit delay

N = length(T) - 1; % don't count unit delay
Ny = 5;
N1 = 1;
Nu = 1;
lambda = 1;

[dmc] = dmc_setup(N, N1, Ny, Nu, lambda, Y11)

k0 = 1+max([1+nb,1+na]);
kmax = tmax/ts + k0;

minLength = length([k0:kmax]) + k0 + Ny;

yr = zeros(minLength, 1);
yr(1:30) = 0.01;
yr(31:end) = -0.05;
yr(61:end) = 0.06;
u = zeros(minLength, 1);
v = zeros(minLength, 1);
vc = zeros(minLength, 1);
y = zeros(minLength, 1);
dvp = zeros(N-N1, 1);

for k = k0:kmax
    % process simulation
    % hammerstein non linearity
    v(k-1) = 1.04*u(k-1) - 14.11*u(k-1)^2 - 16.72*u(k-1)^3 + 562.75*u(k-1)^4;
    
    % ideal scenario
    if IGNORE_NONLINEARITY
        v(k-1) = vc(k-1); % linear bypass (pseudo linear control)
    end
    
    % linear dynamics
    y(k) = B*v(k-1:-1:k-(nb+1)) - A*y(k-1:-1:k-na);
    
    % control action
    dv = dmc_controlaction_constrained(dmc, yr(k+[1:Ny]), dvp, y(k), vc(k-1), [min(vss) max(vss)]);
    
    % pseudo linear control
    vc(k) = vc(k-1) + dv;
    
    % inverse of hammerstein non linearity
    u_roots = flhiInterpolateInverse(interpolant, vc(k));
    
    % iterative search of the root which results in the smallest delta_u
    if(isempty(u_roots))
        warning('DMC HAMMERSTEIN: no solution found');
        u(k) = u(k-1);
    else
        u_roots_diff = abs(u_roots - repmat(u(k-1),size(u_roots,1),1));
        [maxu, maxu_idx] = min(u_roots_diff);
        u(k) = u_roots(maxu_idx);
    end
    
    % shift the vector of past control increments
    dvp(2:end) = dvp(1:end-1);
    dvp(1) = dv;
end

% Integral Squared Error
ISE = sum((yr(k0:kmax)-y(k0:kmax)).^2);
% Integral Squared Variation of Control
ISVC = sum((u(k0+1:kmax)-u(k0:kmax-1)).^2);
J = ISE + lambda*ISVC;
fprintf('ISE = %.4g\n', ISE);
fprintf('ISVC = %.4g\n', ISVC);
fprintf('Cost function J = %.4g\n\n', J);

fontsize = 12;
h = figure('Position', [200 200 400 300]);
subplot(2,1,1);
stairs(([k0:kmax]' - k0)*ts, yr(k0:kmax), 'k--', 'linewidth', 1.5);
hold all;
stairs(([k0:kmax]' - k0)*ts, y(k0:kmax), 'k', 'linewidth', 1.5);
ylabel('Output', 'fontsize', fontsize);
set(gca, 'fontsize', fontsize);
if LIMIT_AXES
    axis(gca, [0 tmax OUTPUT_YMIN OUTPUT_YMAX]);
end
subplot(2,1,2); 
if IGNORE_NONLINEARITY % plot pseudo-linear control signal instead of nonlinear in case of pseudo linear control (ideal scenario)
    u = v;
end
stairs(([k0:kmax]' - k0)*ts, u(k0:kmax), 'k', 'linewidth', 1.5);
ylabel('Control action', 'fontsize', fontsize);
xlabel('Samples', 'fontsize', fontsize);
if LIMIT_AXES
    axis([0 tmax CONTROL_YMIN CONTROL_YMAX]);
end
set(gca, 'fontsize', fontsize);

set(gcf, 'Color', 'w');

export_fig(h, ['DMC_SISO_',FIGURE_NAME,'.pdf'], '-q101', '-p0.01');
