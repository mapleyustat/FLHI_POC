clear;
clc;

s = tf('s');

ts = 0.01;

G1 = 1/(s+1); G1d = c2d(G1, ts);
G2 = 10/(s+1); G2d = c2d(G2, ts);
G3 = 1/(s+10); G3d = c2d(G3, ts);
G4 = 10/(s+10); G4d = c2d(G4, ts);

G = [G1, G2; G3, G4];
Gd = [G1d, G2d; G3d, G4d];

% Compute the steady-state gain matrix
Gss = [sum(G1d.num{1})/sum(G1d.den{1}), sum(G2d.num{1})/sum(G2d.den{1});...
       sum(G3d.num{1})/sum(G3d.den{1}), sum(G4d.num{1})/sum(G4d.den{1})]
   
% Analyze condition and eigen values of Gss
disp('Condition number (Gss)');
cond(Gss)

disp('Eigen-values')
eig(Gss)

if(any(eig(Gss) == 0))
    disp('Ill-conditioned static gain matrix - expect an uncontrollable system');
end

disp('RGA matrix')
RGA = inv(Gss') .* Gss

% Test state and output controllability in state-space form

[A,B,C,D] = ssdata(Gd);

% [A,B,C,D] = ssdata(G);
% Q = C'*C;
% R = eye(2);
% K = lqr(A,B,Q,R)
% Ac = [(A-B*K)];
% Bc = [B];
% Cc = [C];
% Dc = [D];
% states = {'x1' 'x2'};
% inputs = {'r1','r2'};
% outputs = {'y1'; 'y2'};
% sys_cl = ss(Ac,Bc,Cc,Dc,'statename',states,'inputname',inputs,'outputname',outputs);
% t = 0:ts:10;
% r = [ones(size(t)); ones(size(t))]; % Systems might need reference
% % scaling, if a system can't be reference scaled (e.g., two inputs with
% % different scales, it's not controllable
% [y,t,x]=lsim(sys_cl,r,t);
% [AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
% return

n = length(A);
m = size(C,1);
r = size(B,2);

co = [];
oco = [];
for i=1:n
    co = [co A^(i-1)*B];
    oco = [oco C*A^(i-1)*B];
end

co
oco

% REMEMBER: being controllable doesn't mean it can hold a state or output
if rank(oco) == m
    disp('is output controllable!');
else
    disp('not output controllable!');
end

if rank(co) == n
    disp('is state controllable!');
else
    disp('not state controllable!');
end



return


ny = 500;
nu = 1;

model = setmpcsignals(Gd,'MV1',1,'MV2',2);
mpcobj = mpc(model,ts,ny,nu);

Tstop = 30;
Tf = round(Tstop/ts);
r = ones(Tf,2);

sim(mpcobj,Tf,r);