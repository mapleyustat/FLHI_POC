clear all;
clc;

%mapping code from https://www.particleincell.com/2012/quad-interpolation

%create our polygon
px = [0, 1, 1, 0];
py = [0, 0.2, 0.8, 1];
 
%compute coefficients
A=[1 0 0 0;1 1 0 0;1 1 1 1;1 0 1 0];
AI = inv(A);
aa = AI*px';
bb = AI*py';

%x' = (a*l + b*m + c) / (g*l + h*m + 1)       (6)
%y' = (d*l + e*m + f) / (g*l + h*m + 1)
syms a b c d e f g h
ret = solve(0 == c, 0 == f, 1 == (a+c)/(g+1), 0.2 == (d+f)/(g+1), 1 == (a+b+c)/(g+h+1), 0.8 == (d+e+f)/(g+h+1), 0 == (b+c)/(h+1), 1 == (e+f)/(h+1))

% CBimage = checkerboard(50) > 0.5;
CBimage = zeros(400);
for i=0:399
    for j=0:399
        CBimage(i+1,j+1) = (i/399) * (j/399);
    end
end


max_width = size(CBimage, 2);
max_height = size(CBimage, 1);

CBinterp = ones(max_height, max_width);

a = double(vpa(ret.a));
b = double(vpa(ret.b));
c = double(vpa(ret.c));
d = double(vpa(ret.d));
e = double(vpa(ret.e));
f = double(vpa(ret.f));
g = double(vpa(ret.g));
h = double(vpa(ret.h));

for i=0:max_height-1
    for j=0:max_width-1
        l = j/(max_width-1);
        m = i/(max_height-1);
        lm = [1; l; m; l*m];
        %bilinear transform
        x = aa'*lm;
        y = bb'*lm;
        %projective transform
        % this one presents gaps because we're doing the opposite (going
        % from square box to skewed box), so it misses some spots
%         x = (a*l + b*m + c) / (g*l + h*m + 1);
%         y = (d*l + e*m + f) / (g*l + h*m + 1);
        
        x = round(x*(max_width - 1) + 1);
        y = round(y*(max_height - 1) + 1);
        
        CBinterp(y,x) = CBimage(i+1,j+1);
    end
end

figure;
imshow(CBimage);

figure;
imshow(CBinterp)
