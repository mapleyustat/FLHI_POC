clearvars -except FIGURE_NAME INTERPOLANT_TYPE OUTPUT_YMAX OUTPUT_YMIN CONTROL_YMAX CONTROL_YMIN LIMIT_AXES IGNORE_NONLINEARITY

disp('====================');
fprintf('=  DMC MIMO %s\n', FIGURE_NAME);
disp('====================');

ts = 1;
tmax = 120;
N = 31;
T = 0:ts:(N-1)*ts;

Gd11 = filt([0, 0.1, 0.2], [1, -1.2, 0.35], ts);
Gd12 = filt([0, 1], [1, -0.7], ts);
Gd21 = filt([0, 0.3, 0.2], [1, -0.8], ts);
Gd22 = filt([0, 1, 0.5], [1, 0.4], ts);

Gss = [sum(Gd11.num{1})/sum(Gd11.den{1}), sum(Gd12.num{1})/sum(Gd12.den{1}); sum(Gd21.num{1})/sum(Gd21.den{1}), sum(Gd22.num{1})/sum(Gd22.den{1})];
RGA = (Gss.')^-1.*Gss

uss1 = [-2:0.2:2];
uss2 = [-2:0.2:2];
vss1 = uss1.^3 - uss1.*uss2 + 2*uss2.^2;
vss2 = 0.582*(exp(uss1 + uss2) - 1);

% [X,Y] = meshgrid(uss1, uss2);
% Z1 = X.^3 - X.*Y + 2*Y.^2;
% Z2 = 0.582*(exp(X+Y) - 1);
% figure;
% surf(X,Y,Z1);
% figure;
% surf(X,Y,Z2);

pointVector = [];

for i = 1:length(uss1)
    for j = 1:length(uss2)
        point.xi = [uss1(i), uss2(j)];
        point.xo = [vss1(i), vss2(j)];
    
        pointVector = [pointVector; point];
    end
end

[interpolant] = flhiInterpolant(pointVector, INTERPOLANT_TYPE)

B11 = conv(Gd11.num{1}(2:end), Gd12.den{1});
B12 = conv(Gd12.num{1}(2:end), Gd11.den{1});
A1 = conv(Gd11.den{1}, Gd12.den{1});
A1 = A1(2:end);

B21 = conv(Gd21.num{1}(2:end), Gd22.den{1});
B22 = conv(Gd22.num{1}(2:end), Gd21.den{1});
A2 = conv(Gd21.den{1}, Gd22.den{1});
A2 = A2(2:end);

nb1 = length(B11) - 1;
nb2 = length(B21) - 1;
na1 = length(A1);
na2 = length(A2);

Y11 = step(Gd11, T); 
Y12 = step(Gd12, T); 
Y21 = step(Gd21, T); 
Y22 = step(Gd22, T); 

input_count = 2;
output_count = 2;
Y = zeros(length(T) - 1, input_count, output_count); % skip unit delay
Y(:,:,1) = [Y11(2:end), Y12(2:end)];
Y(:,:,2) = [Y21(2:end), Y22(2:end)];
 
N = length(T) - 1; % don't count unit delay
Ny = 15;
N1 = 1;
Nu = 5;
lambda = 5;

[dmc] = dmc_setup(N, N1, Ny, Nu, lambda, Y)

k0 = 1+max([1+nb1,1+na1,1+nb2,1+na2]);
kmax = tmax/ts + k0;

minLength = length([k0:kmax]) + k0 + Ny;

yr = 0.1*ones(minLength, output_count);
yr(1:40, 1) = 11;
yr(41:end, 1) = 15;
yr(1:40, 2) = 6;
yr(41:80, 2) = 10;
yr(81:end, 2) = 13;
yr(111:end, 2) = 10;
u = zeros(minLength, input_count);
v = zeros(minLength, input_count);
vc = zeros(minLength, input_count);
y = zeros(minLength, output_count);
f = zeros(minLength, output_count);
dvp = zeros(N-N1, input_count);


for k = k0:kmax
    % process simulation
    % hammerstein non linearity
    v(k-1, 1) = u(k-1, 1).^3 - u(k-1, 1).*u(k-1, 2) + 2*u(k-1, 2).^2;
    v(k-1, 2) = 0.582*(exp(u(k-1, 1)+u(k-1, 2)) - 1);
    
    % ideal scenario
    if IGNORE_NONLINEARITY
        v(k-1,:) = vc(k-1,:); % linear bypass (pseudo linear control)
    end
    
    % linear dynamics
    y(k, 1) = B11*v(k-1:-1:k-(nb1+1), 1) + B12*v(k-1:-1:k-(nb1+1), 2) - A1*y(k-1:-1:k-na1, 1);
    y(k, 2) = B21*v(k-1:-1:k-(nb2+1), 1) + B22*v(k-1:-1:k-(nb2+1), 2) - A2*y(k-1:-1:k-na2, 2);
    
    % control action
    dv = dmc_controlaction_constrained(dmc, yr(k+[1:Ny-N1+1], :), dvp, y(k,:), vc(k-1, :), [min(vss1) max(vss1); min(vss2) max(vss2)]);
    %dv = dmc_controlaction(dmc, yr(k+[1:Ny-N1+1], :), dvp, y(k,:));
    
    % pseudo linear control
    vc(k, 1) = vc(k-1, 1) + dv(1);
    vc(k, 2) = vc(k-1, 2) + dv(2);
    
    % inverse of hammerstein non linearity
    u_roots = flhiInterpolateInverse(interpolant, vc(k, :));
    
    % iterative search of the root which results in the smallest delta_u
    if(isempty(u_roots))
        warning('DMC HAMMERSTEIN: no solution found');
        u(k, :) = u(k-1, :);
    else
        u_roots_diff = abs(u_roots - repmat(u(k-1, :),size(u_roots,1),1));
        [maxu, maxu_idx] = min(sum(u_roots_diff,2));
        u(k, :) = u_roots(maxu_idx, :);
    end
    
    % shift the vector of past control increments
    dvp(2:end,:) = dvp(1:end-1,:);
    dvp(1,:) = dv';
end

% figure;
% subplot(2,1,1);
% stairs(([k0:kmax]' - k0)*ts, yr(k0:kmax,:), 'k--');
% hold all;
% ylabel('Sa�da (mag)');
% xlabel('Tempo (t)');
% stairs(([k0:kmax]' - k0)*ts, y(k0:kmax,:), 'LineWidth', 2);
% subplot(2,1,2);
% stairs(([k0:kmax]' - k0)*ts, u(k0:kmax, :), 'LineWidth', 2);
% ylabel('Controle (mag)');
% xlabel('Tempo (t)');
% figure;
% subplot(3,1,1);
% stairs(([k0:kmax]' - k0)*ts, v(k0:kmax,:), 'LineWidth', 2);
% subplot(3,1,2);
% stairs(([k0:kmax]' - k0)*ts, vc(k0:kmax, :), 'LineWidth', 2);
% subplot(3,1,3);
% stairs(([k0:kmax]' - k0)*ts, abs(v(k0:kmax, :) - vc(k0:kmax, :)), 'LineWidth', 2);

ISE = sum((yr(k0:kmax,1)-y(k0:kmax,1)).^2) + sum((yr(k0:kmax,2)-y(k0:kmax,2)).^2);
ISVC = sum((u(k0+1:kmax,1)-u(k0:kmax-1,1)).^2) + sum((u(k0+1:kmax,2)-u(k0:kmax-1,2)).^2);
J = ISE + lambda*ISVC;
fprintf('ISE = %.4g\n', ISE);
fprintf('ISVC = %.4g\n', ISVC);
fprintf('Cost function J = %.4g\n\n', J);

fontsize = 12;
h = figure('Position', [200 200 400 300]);
subplot(2,1,1);
stairs(([k0:kmax]' - k0)*ts, yr(k0:kmax,:), 'k--', 'linewidth', 1.5);
hold all;
stairs(([k0:kmax]' - k0)*ts, y(k0:kmax,1), 'k', 'linewidth', 1.5);
hold all;
stairs(([k0:kmax]' - k0)*ts, y(k0:kmax,2), 'linewidth', 1.5, 'color', [0,0,0]+0.5);
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
stairs(([k0:kmax]' - k0)*ts, u(k0:kmax,2), 'linewidth', 1.5, 'color', [0,0,0]+0.5);
ylabel('Control action', 'fontsize', fontsize);
xlabel('Samples', 'fontsize', fontsize);
if LIMIT_AXES
    axis([0 tmax CONTROL_YMIN CONTROL_YMAX]);
end
set(gca, 'fontsize', fontsize);

set(gcf, 'Color', 'w');

export_fig(h, ['DMC_MIMO_',FIGURE_NAME,'.pdf'], '-q101', '-p0.01');


