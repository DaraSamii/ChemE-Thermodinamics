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
t = [300:0.1:310];

BUBL_P = zeros(size(X, 2), size(t, 2));
DEW_P = zeros(size(X, 2), size(t, 2));
for i = 1:1:size(X,2)
    for j = 1:1:size(t, 2)
        T = t(j);
        x = [X(i) 1-X(i)];
        y = x;
        
        BUBL_P(i,j) = bubble_P_gammaphi(x, T, Tc, w, Zc, Vc, Pc, a, b, c, par, a_mn);
        DEW_P(i,j) = dew_P_gammaphi(y, T, Tc, w, Zc, Vc, Pc, a, b, c, par, a_mn);
    end
end
%% Ploting 

[X,T] = meshgrid(X, t);

figure(1);
surf(X, T, DEW_P');
hold on;
mesh(X, T, BUBL_P');

title('Pressure acording to x,y and Temperature');
xlabel('x,y');
ylabel('Temperature');
zlabel('Pressure');
