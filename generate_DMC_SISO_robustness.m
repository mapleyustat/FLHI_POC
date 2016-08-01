clearvars -except FIGURE_NAME INTERPOLANT_TYPE

disp('====================');
fprintf('=  DMC SISO ROBUSTNESS %s\n', FIGURE_NAME);
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

% % % % % % % % %
% DMC SISO control action without input prediction is given by:
% d*u = K*(S*yr - dG*d*u - S*y)
% where S = ones(Ny,1).
%
% Closed loop transfer function can be obtained by the following code:
% clear all; clc;
% syms S u yr y B A X K d dG
% eq1 = d*u == K*(S*yr - dG*d*u - S*y)
% eq2 = y == (B/A)*u;
% eq3 = y == X*yr;
% res = solve([eq1,eq2,eq3], [u,y,X])
%
% such that res.X is the closed loop transfer function:
% res.X = (B*K*S)/(A*d + B*K*S + A*K*dG*d)
%
% from the denominator it's easy to see the open loop transfer function L is:
% L = (B*K*S)/(A*d + A*K*dG*d)
% % % % % % % % %

% Gain margin, phase margin and sensitivity calculation from above
B = filt([0 B], 1, ts);
A = filt([1 A], 1, ts);
d = filt([1 -1], 1, ts);
S = ones(dmc.Ny, 1);
KS = filt(dmc.K*S, 1, ts); % sum of K
KdG = filt([0 dmc.K*dmc.dG], 1, ts); % dG is terms of PAST u inputs
%CL = B*KS/(A*d + B*KS + A*KdG*d);
L = B*KS/(A*d + A*KdG*d);

% % this code is equivalent to the above but in RST terms
% R=filt([1 dmc.K*dmc.dG],1,ts)*d;
% T=filt(sum(dmc.K),1,ts);
% S=filt(sum(dmc.K),1,ts);
% CL_equivalent = B*T/(A*R+B*S)
% L_equivalent = B*S/(A*R)
% L_equivalent2 = (1/R)*(B/A)*S

disp('Robustness information:');
[Gm,Pm,Wgm,Wpm] = margin(L);
Gm
Pm

% This should be equal to Ms = norm(1/(1+L),inf)
wcst = wcsens(L);
Ms = wcst.Si.MaximumGain.UpperBound

%% Maximize the absolute error between FLHI and true nonlinearity

disp('Maximum absolute error between models:');

% true nonlinearity
NL = @(uss)1.04.*uss - 14.11.*uss.^2 - 16.72.*uss.^3 + 562.75.*uss.^4;

% flhi
FLHI = @(xi) flhiInterpolate(interpolant, xi);

% objective function - maximization of model absolute error (MAE) or model
% absolute relative error (MARE)
modelError = @(xi) -abs(NL(FLHI(xi))/xi); % MARE
%modelError = @(xi) -abs(NL(FLHI(xi)) - xi); % MAE

% optimize the above problem in terms of minimization
%optoptions = optimset('Display','off','Algorithm','interior-point','TolFun',1e-10);
%[xi,fval] = fmincon(modelError,[0],[],[],[],[],[-0.15],[0.15],[],optoptions);

optoptions = optimset('Display','off','Algorithm','interior-point','TolFun',1e-10);
x0 = [0.001]; % start point away from the minimum
lb = [-0.15];
ub = [0.15];
problem = createOptimProblem('fmincon','objective',modelError,'x0',x0,'lb',lb,'ub',ub,'options',optoptions);
%gs = GlobalSearch;
%[xg,fg,flg,og] = run(gs,problem)
ms = MultiStart;
[xg,fg,flg,og] = run(ms,problem,100);

% turn minimization into maximization
xg
fg = -fg
