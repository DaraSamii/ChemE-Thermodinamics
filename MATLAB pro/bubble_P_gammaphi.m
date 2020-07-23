function [b_P,y] = bubble_P_gammaphi(x, T, Tc, w, Zc, Vc, Pc, a, b, c, par, a_mn)
%bubble_P_gammaphi calculates the bubble pressure and vapor phase
%composition of a liquid mixture
 
%INPUT:
% x - row vector of the mole fraction of the components
% T - Temperature in Kelvin
% Tc - row vector of the critical temperature (K) of each species
% w - row vector 
% Zc - row vector
% Vc - row vector of the critical volume (cm3/mol) of each species
% Pc - row vector of the critical pressure (bar) of each species
% a,b,c - row vectors containing the Antoine constants of each species
% par, a_mn - parameters for UNIFAC
 
%OUTPUT
% y - vapor phase composition of the liquid mixture
% b_P - bubble pressure in bar
 
%matrix holders for the values of PHI, Psat, z, and y
PHI = ones (size(x,2),1);
Psat = zeros (size(x,2),1);
z = zeros (1,size(x,2));
y = zeros (1,size(x,2));
%solves the gamma of the species using UNIFAC
gamma = gamma_UNIFAC(x, T, par, a_mn);
%calculates the saturated pressure (kPa) of each species for a given temperature
for i = 1:size(x,2)
    Psat(i) = exp(a(i)-(b(i)/((T-273.15)+c(i))));
end
tol = 10^-6; %tolerance for the iterations
diff = 1;
%iterations for the calculation of bubble P and vapor composition
fprintf('\n\nIteration values of P (bubble P) in bar :\n\n');
while abs(diff)>tol
    %for loop that calculates the value of P
    for i = 1:size(x,2)
        z(i) = (x(i)*gamma(i)*Psat(i))/PHI(i);
    end
    P = sum(z);
    %calculates the value of Bij
    B = B_ij(T, Tc, w, Zc, Vc, Pc);
    %for loop that calculates the vapor phase composition of each species
    for i = 1:size(x,2)
        y(i) = (x(i)*gamma(i)*Psat(i))/(P*PHI(i));
    end
    %calculates phihat and phisat for each species
    phi = phi_i(y, B, T, P, Psat);
    %for loop that calculates the value of PHI and the bubble pressure
    for i = 1:size(x,2)
        PHI(i) = (phi(i,1)/phi(i,2));
        z(i) = (x(i)*gamma(i)*Psat(i))/PHI(i);
    end
    Pnew = sum(z); %bubble pressure in kPa
    diff = P-Pnew;
    disp(P/100);
end
b_P = Pnew/100; %bubble pressure in bar
end
