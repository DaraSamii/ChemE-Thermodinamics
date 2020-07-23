close all;
clear;
clc;

%% Constant Variables
Tc = [506.6 512.6];
w = [0.331 0.564];
Zc = [0.257 0.224];
Vc = [228 118];
Pc = [47.5 80.97];
a = [14.2456 16.5785];
b = [2662.78 3638.27];
c = [219.69 239.5];
par = [1 1 0.9011 0.848 1 0;11 21 1.9031 1.728 1 0;6 15 1.4311 1.432 0 1];
a_mn = [0 232.1 697.2;114.8 0 249.6;16.51 -10.72 0];

%% calculating DEW Pressure and Bubble Pressure 

X = [0:.1:1];
p = [300:0.1:310];

BUBL_T = zeros(size(X, 2), size(p, 2));
DEW_T = zeros(size(X, 2), size(p, 2));

for i = 1:1:size(X,2)
    for j = 1:1:size(p, 2)
        P = p(j);
        x = [X(i) 1-X(i)];
        y = x;
        
        BUBL_T(i, j) = bubble_T_gammaphi(x, Tc, w, Zc, Vc, P, Pc, a, b, c, par, a_mn);
        DEW_T(i, j) = dew_T_gammaphi(y, P, Tc, w, Zc, Vc, Pc, a, b, c, par, a_mn);
    end
end
%% Ploting 

[X, P] = meshgrid(X, p);

figure(1);
surf(X, P, DEW_T');
hold on;
mesh(X, P, BUBL_T');

title('Temperature acording to x,y and Pressure');
xlabel('x,y');
ylabel('Pressure');
zlabel('Temperature');
