function [d_P,x] = dew_P_gammaphi(y,T,Tc,w,Zc,Vc,Pc,a,b,c,par,a_mn)
%dew_P_gammaphi calculates the dew pressure and liquid phase
%composition of a liquid mixture
 
%INPUT:
%y - row vector of the mole fraction of the components
%T - Temperature in Kelvin
%Tc - row vector of the critical temperature (K) of each species
%w - row vector
%Zc - row vector
%Vc - row vector of the critical volume (cm3/mol) of each species
%Pc - row vector of the critical pressure (bar) of each species
%a,b,c - row vectors containing the Antoine constants of each species
%par, a_mn - parameters for UNIFAC
 
%OUTPUT
%x - liquid phase composition of the mixture
%d_P - dew pressure in bar
 
%matrix holders for the values of PHI, gamma, Psat, Z, X, and Xnew
PHI = ones(size(y,2),1); %PHI = 1
gamma =[1,1]; %gamma = 1
Psat = zeros(size(y,2),1);
z = zeros(1,size(y,2));
x = zeros(1,size(y,2));
xnew = zeros(1,size(y,2));
tol = 10^-6; %tolerance
%solves for the saturated pressure (kPa) of each species and the Pressure
%where PHI = 1 and gamma = 1
for i = 1:size(y,2)
    Psat(i) = exp(a(i)-(b(i)/((T-273.15)+c(i))));
    z(i) = (y(i)*PHI(i))/(gamma(i)*Psat(i));
end
P = 1/sum(z);
%calculates the new liquid phase composition using the pressure obtained
for i = 1:size(y,2)
    x(i) = (y(i)*P*PHI(i))/(gamma(i)*Psat(i));
end
%solves for the value of gamma for each species
%bubble_P_gammaphi 
%for loop that calculates the new pressure using the new values of gamma
%and PHI = 1;
for i = 1:size(y,2)
    z(i) = (x(i)*gamma(i)*Psat(i))/PHI(i);
end
P = sum(z);
 
error = 1;
fprintf('\n\nIteration values of P (dew P) in kPa :\n\n');
while abs(error) > tol
    B = B_ij(T, Tc, w, Zc, Vc, Pc); %solves for the values of Bij
    phi = phi_i(y, B, T, P, Psat); %solves for phihat and phisat
    diff = 1;
    while abs(diff) > tol
        %calculates the new value of PHI and x
        for i = 1:size(y,2)
            PHI(i) = (phi(i,1)/phi(i,2));
            x(i) = (y(i)*PHI(i)*P)/(gamma(i)*Psat(i));
        end
        %normalizes the value of x
        while sum(x) ~= 1
            for i = 1:size(y,2)
                xnew(i) = x(i)/sum(x);
            end
            x = xnew;
        end
        %calculates the new value of gamma using the new value of x
        gammanew = gamma_UNIFAC(x, T, par, a_mn);
        disp(gammanew);
        diff = gammanew - gamma;
        gamma = gammanew;
    end
    %calculates the new Pressure using the new values of PHI, gamma and x
    for i = 1:size(y,2)
        z(i) = (x(i)*gamma(i)*Psat(i))/PHI(i);
    end
    Pnew = sum(z);
    error = Pnew-P;
    P = Pnew;
    disp(P);
end
d_P = P/100;
end
