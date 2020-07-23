close all;
clear;
clc;


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

%% Calculating BUBL_P and y

x = [0.3 0.7];
T = 348.15;

[BUBL_P, B_P_y] = bubble_P_gammaphi(x, T, Tc, w, Zc, Vc, Pc, a, b, c, par, a_mn);

%% Calculating DEW_P and x

y = [0.4309 0.591];
T = 348.15;

[DEW_P, D_P_x] = dew_P_gammaphi(y, T, Tc, w, Zc, Vc, Pc, a, b, c, par, a_mn);

%% Calculating BUBL_T and y

x = [0.3 0.7];
P = 2.0144;

[BUBL_T, B_T_y] = bubble_T_gammaphi(x, Tc, w, Zc, Vc, P, Pc, a, b, c, par, a_mn);

%% Calcuating DEW_T and x

y = [0.4309 0.5691];
P = 2.0144;

[DEW_T , D_T_x] = dew_T_gammaphi(y, P, Tc, w, Zc, Vc, Pc, a, b, c, par, a_mn);

%% Displaying Results
disp(['at T = 348.15 and x1 = 0.3 and x2 = 0.7 Bubble Pressure is ',...
        num2str(BUBL_P), ' and y1 = ',num2str(B_P_y(1)), ' and y2 = ' , num2str(B_P_y(2)),'.']);

disp(['at P = 2.0144 and x1 = 0.3 and x2 = 0.7 Buble Temperature is ',...
        num2str(BUBL_T), ' and y1 = ',num2str(B_T_y(1)), ' and y2 = ' , num2str(B_T_y(2)),'.']);
    
disp(['at T = 348.15 and y1 = 0.4309 and y2 = 0.5691 Dew Pressure is ',...
        num2str(DEW_P), ' and x1 = ',num2str(D_P_x(1)), ' and x2 = ' , num2str(D_P_x(2)),'.']);
    
disp(['at P = 2.0144 and y1 = 0.4309 and y2 = 0.5691 Dew Temperature is ',...
        num2str(DEW_T), ' and x1 = ',num2str(D_T_x(1)), ' and x2 = ' , num2str(D_T_x(2)),'.']);









