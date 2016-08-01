clear all;
clc;

%mapping code from https://www.particleincell.com/2012/quad-interpolation

%create our polygon
px = [0, 1, 1, 0];
py = [0, 0.1, 0.9, 1];
 
%compute coefficients
A=[1 0 0 0;1 1 0 0;1 1 1 1;1 0 1 0];
AI = inv(A);
a = AI*px';
b = AI*py';

CBimage = checkerboard(50) > 0.5;

max_width = size(CBimage, 2);
max_height = size(CBimage, 1);

CBinterp = ones(max_height, max_width);

for i=1:max_height
    for j=1:max_width
        l = j/max_width;
        m = i/max_height;
        lm = [1; l; m; l*m];
        x = round((a'*lm)*(max_width));
        y = round((b'*lm)*(max_height));
        
        CBinterp(y,x) = CBimage(i,j);
    end
end

figure;
imshow(CBimage);

figure;
imshow(CBinterp)
