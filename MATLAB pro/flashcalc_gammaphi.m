function [x,y] = flashcalc_gammaphi(T,P,z,Tc,w,Zc,Vc,Pc,a,b,c,par,a_mn)
%flashcalc_gammaphi calculates the vapor phase and liquid phase
%compositions of a liquid mixture
 
%INPUT:
%P - Pressure in kPa
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
%y - vapor phase composition of the mixture
 
PG = P; 
%matrix holders for the values of xold, yold, Psat, PHI, K, u, q and gamma
xold = zeros(1,size(z,2));
yold = zeros(1,size(z,2));
Psat = zeros(size(z,2),1);
PHI = zeros(size(z,2),1);
K = zeros(size(z,2),1);
u = zeros(1,size(z,2));
q = zeros(1,size(z,2));
gamma = zeros(size(z,2),1);
%calculates the saturated pressure (kPa) of each species
for i = 1:size(z,2)
    Psat(i) = exp(a(i)-(b(i)/((T-273.15)+c(i))));
end
%for the dew point values
y = z; %vapor phase composition
V_d = 1; %value of V for dew point
%calculates the dew pressure (bar) and liquid phase composition
[d_P,x] = dew_P_gammaphi(y,T,Tc,w,Zc,Vc,Pc,a,b,c,par,a_mn); 
gamma_d = gamma_UNIFAC(x, T, par, a_mn); %calculates the gamma of each species
B = B_ij(T, Tc, w, Zc, Vc, Pc); %calculates the Bij 
P = d_P*100; %converts the dew pressure to kPa
phi = phi_i(y, B, T, P, Psat); %give the values of phihat and phisat
%calculates the values of PHI
for i = 1:size(z,2)
    PHI(i) = (phi(i,1)/phi(i,2));
end
PHI_d = PHI;
 
%for bubble point values
x = z; %liquid phase composition
V_b = 0; %value of V for bubble point
%calculates the bubble pressure (bar) and vapor phase composition
[b_P,y] = bubble_P_gammaphi(x,T,Tc,w,Zc,Vc,Pc,a,b,c,par,a_mn);
gamma_b = gamma_UNIFAC(x, T, par, a_mn); %calculates the gamma of each species
B = B_ij(T, Tc, w, Zc, Vc, Pc); %calculates  the Bij
P = b_P*100; %converts the bubble pressure to kPa
phi = phi_i(y, B, T, P, Psat); %give the values of phihat and phisat
%calculates the values of PHI
for i = 1:size(z,2)
    PHI(i) = (phi(i,1)/phi(i,2));
end
PHI_b = PHI;
 
%for flash calculation values
PG = PG/100; %converts the given pressure to bar
if PG < b_P && PG > d_P 
    %interpolates the values of V, gamma and PHI from the calculates values
    %of V, gamma, and PHI from the dew and bubble point values
    V = V_d - ((d_P-P)*(V_d-V_b))/(d_P-b_P);
    for i = 1:size(z,2)
        gamma(i) = gamma_d(i) - ((d_P-P)*(gamma_d(i)-gamma_b(i)))/(d_P-b_P);
        PHI(i) = PHI_d(i) - ((d_P)*(PHI_d(i)-PHI_b(i))/(d_P-b_P));
    end
    diff1 = 1;
    diff2 = 1;
    diff3 = 1;
    tol = 10^-6;
    PG = PG*100; %converts the given pressure to kPa
    while abs(diff1) > tol && abs(diff2) > tol && abs(diff3) > tol
        %calculates the values of K, F (sum of u) and dF/dV (sum of q)
        for i = 1:size(z,2)
            K(i) = (gamma(i)*Psat(i))/(PHI(i)*PG);
            u(i) = (z(i)*(K(i)-1))/(1+V*(K(i)-1));
            q(i) = (z(i)*(K(i)-1)^2)/((1+V*(K(i)-1))^2);
        end
        F = sum(u);
        dFdV = -sum(q);
        Vnew = V - (F/dFdV); %new value of V using Newton's method
        %solves for the vapor and liquid phase compositions
        for i = 1:size(z,2)
            x(i) = z(i)/(1+Vnew*(K(i)-1));
            y(i) = K(i)*x(i);
        end
        gamma = gamma_UNIFAC(x, T, par, a_mn); %solves for the gamma of each species 
        B = B_ij(T, Tc, w, Zc, Vc, Pc); %solves for Bij
        P = PG;
        phi = phi_i(y, B, T, P, Psat); %solves for phihat and phisat
        %solves for PHI
        for i = 1:size(z,2)
            PHI(i) = (phi(i,1)/phi(i,2));
        end
        PG = P;
        diff1 = V - Vnew;
        V = Vnew;
        diff2 = xold - x;
        xold = x;
        diff3 = yold - y;
        yold = y;
    end
    x = xold;
    y = yold;
end
end
