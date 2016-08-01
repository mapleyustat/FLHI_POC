clearvars -except FIGURE_NAME INTERPOLANT_TYPE OUTPUT_YMAX OUTPUT_YMIN CONTROL_YMAX CONTROL_YMIN LIMIT_AXES IGNORE_NONLINEARITY

disp('====================');
fprintf('=  DMC MISO UNCOUPLED %s\n', FIGURE_NAME);
disp('====================');

ts = 1;
tmax = 100;
N = 31;
T = 0:ts:(N-1)*ts;

Gd1 = filt([0, 0.1333, 0.0667], [1, -1.5, 0.7], ts);
Gd2 = filt([0, 0.4, 0.3], [1, -0.9, 0.6], ts);

Gss = [sum(Gd1.num{1})/sum(Gd1.den{1}), sum(Gd2.num{1})/sum(Gd2.den{1})];
RGA = pinv(Gss.').*Gss

uss1 = [-2:0.1:2];
uss2 = [-2:0.1:2];
vss1 = uss1 + 4*uss1.^2 + 1.5*uss1.^3;
vss2 = uss2 + 3*uss2.^2 + 2*uss2.^3;

pointVector1 = [];
pointVector2 = [];

for i = 1:length(uss1)
    point.xi = uss1(i);
    point.xo = vss1(i);
    
    pointVector1 = [pointVector1; point];
    
    point.xi = uss2(i);
    point.xo = vss2(i);
    
    pointVector2 = [pointVector2; point];
end

[interpolant1] = flhiInterpolant(pointVector1, INTERPOLANT_TYPE)
[interpolant2] = flhiInterpolant(pointVector2, INTERPOLANT_TYPE);

B1 = conv(Gd1.num{1}(2:end), Gd2.den{1});
B2 = conv(Gd2.num{1}(2:end), Gd1.den{1});
A = conv(Gd1.den{1}, Gd2.den{1});
A = A(2:end);

nb = length(B1) - 1;
na = length(A);

Y11 = step(Gd1, T); 
Y12 = step(Gd2, T); 

% skip unit delay
Y = [Y11(2:end), Y12(2:end)];

N = length(T) - 1; % don't count unit delay
Ny = 10;
N1 = 1;
Nu = 5;
lambda = 1;

[dmc] = dmc_setup(N, N1, Ny, Nu, lambda, Y)

k0 = 1+max([1+nb,1+na]);
kmax = tmax/ts + k0;

minLength = length([k0:kmax]) + k0 + Ny;

yr = 5*ones(minLength, 1);
yr(1:30) = 10;
yr(31:end) = 50;
yr(61:end) = -5;
u = zeros(minLength, 2);
v = zeros(minLength, 2);
vc = zeros(minLength, 2);
y = zeros(minLength, 1);
dvp = zeros(N-N1, 2);

for k = k0:kmax
    % process simulation
    % hammerstein non linearity
    v(k-1, 1) = u(k-1, 1) + 4*u(k-1, 1)^2 + 1.5*u(k-1, 1)^3;
    v(k-1, 2) = u(k-1, 2) + 3*u(k-1, 2)^2 + 2*u(k-1, 2)^3;
    
    % ideal scenario
    if IGNORE_NONLINEARITY
        v(k-1,:) = vc(k-1,:); % linear bypass (pseudo linear control)
    end
    
    % linear dynamics
    y(k) = B1*v(k-1:-1:k-(nb+1), 1) + B2*v(k-1:-1:k-(nb+1), 2) - A*y(k-1:-1:k-na);
    
    % control action
    dv = dmc_controlaction_constrained(dmc, yr(k+[1:Ny]), dvp, y(k), vc(k-1, :), [min(vss1) max(vss1); min(vss2) max(vss2)]);
    
    % pseudo linear control
    vc(k, 1) = vc(k-1, 1) + dv(1);
    vc(k, 2) = vc(k-1, 2) + dv(2);
    
    % inverse of hammerstein non linearity
    u_roots = flhiInterpolateInverse(interpolant1, vc(k, 1));
    
    % iterative search of the root which results in the smallest delta_u
    if(isempty(u_roots))
        warning('DMC HAMMERSTEIN 1: no solution found');
        u(k, 1) = u(k-1, 1);
    else
        u_roots_diff = abs(u_roots - repmat(u(k-1, 1),size(u_roots,1),1));
        [maxu, maxu_idx] = min(u_roots_diff);
        u(k, 1) = u_roots(maxu_idx);
    end
    
    u_roots = flhiInterpolateInverse(interpolant2, vc(k, 2));
    if(isempty(u_roots))
        warning('DMC HAMMERSTEIN 2: no solution found');
        u(k, 2) = u(k-1, 2);
    else
        u_roots_diff = abs(u_roots - repmat(u(k-1, 1),size(u_roots,1),1));
        [maxu, maxu_idx] = min(u_roots_diff);
        u(k, 2) = u_roots(maxu_idx);
    end
    
    % shift the vector of past control increments
    dvp(2:end,:) = dvp(1:end-1,:);
    dvp(1,:) = dv';
end

ISE = sum((yr(k0:kmax)-y(k0:kmax)).^2);
ISVC = sum((u(k0+1:kmax,1)-u(k0:kmax-1,1)).^2) + sum((u(k0+1:kmax,2)-u(k0:kmax-1,2)).^2);
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
stairs(([k0:kmax]' - k0)*ts, u(k0:kmax,1), 'k', 'linewidth', 1.5);
hold all;
stairs(([k0:kmax]' - k0)*ts, u(k0:kmax,2), 'k--', 'linewidth', 1.5);
ylabel('Control action', 'fontsize', fontsize);
xlabel('Samples', 'fontsize', fontsize);
if LIMIT_AXES
    axis([0 tmax CONTROL_YMIN CONTROL_YMAX]);
end
set(gca, 'fontsize', fontsize);

set(gcf, 'Color', 'w');

export_fig(h, ['DMC_MISO_UNCOUPLED_',FIGURE_NAME,'.pdf'], '-q101', '-p0.01');
