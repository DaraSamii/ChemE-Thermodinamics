function [d_T,x] = dew_T_gammaphi(y,P,Tc,w,Zc,Vc,Pc,a,b,c,par,a_mn)
%dew_T_gammaphi calculates for the composition of the mixture in the liquid
%phase and the dew point temperature in K
 
%INPUT:
%y - vapor phase composition of the mixture
%P - Pressure in bar
%Tc - row vector of the critical temperature (K) of each species
%w - row vector 
%Zc - row vector
%Vc - row vector of the critical volume (cm3/mol) of each species
%Pc - row vector of the critical pressure (bar) of each species
%A,b,C - row vectors containing the Antoine constants of each species
%par, a_mn - parameters for UNIFAC
 
%OUTPUT
%x - liquid phase composition of the mixture
%T -  dew point temperature in K
 
PHI = ones(size(y,2),1);
gamma = ones(size(y,2),1);
Tsat = zeros(size(y,2),1);
Psat = zeros(size(y,2),1);
z = zeros(1,size(y,2));
x = zeros(1,size(y,2));
xnew = zeros(1,size(y,2));
%convert input P from bar to kPa
P = P*100; %P(kPa)
 
%find value of Tsat using the given P
for i=1:size(y,2)
    Tsat(i)=((b(i)/(a(i)-log(P)))-c(i)); %T(oC)
end
 
%find the initial value of T
for i=1:size(y,2)
    z(i)=((y(i)*Tsat(i)));
end
T=sum(z); %T(oC)
 
%find the value of Psat(i) using Reverse Antoine
for i=1:size(y,2)
    Psat(i)=exp(a(i)-(b(i)/(T+c(i)))); %P(kPa)
end
 
%calculating Pjsat
Pjsat=Psat(1);
for i=1:size(y,2) %1x2
    z(i)=((y(i)*PHI(i))/gamma(i))*(Pjsat/Psat(i));
end
j=sum(z);
Pjsat=P*j; %(kPa)
   
%find the new value of T from Pjsat
T=((b(1)/(a(1)-log(Pjsat)))-c(1)); % T(oC)
 
%solve for Psat(i)
for i=1:size(y,2)
    Psat(i)=exp(a(i)-(b(i)/(T+c(i)))); %kPa
end
 
%solve for phi(i) and x
T=T+273.15; %K
B = B_ij(T, Tc, w, Zc, Vc, Pc);
phi = phi_i(y, B, T, P, Psat);
for i = 1:size(y,2)
    PHI(i) = (phi(i,1)/phi(i,2));
    x(i) = (y(i)*PHI(i)*P)/(gamma(i)*Psat(i));
end
 
%solve for gamma
gamma = gamma_UNIFAC(x, T, par, a_mn);
 
%solve for Pjsat
Pjsat=Psat(1);
for i=1:size(y,2) %1x2
    z(i)=((y(i)*PHI(i))/gamma(i))*(Pjsat/Psat(i));
end
j=sum(z);
Pjsat=P*j; %(kPa)
 
%solve for T
T=((b(1)/(a(1)-log(Pjsat)))-c(1)); % T(oC)
 
Tdiff = 1;
tol = 10^-6;
int = 1;
fprintf('\n\nIteration values of T (dew T) in K :\n\n');
while abs(Tdiff)>tol
    %find Psat
    for i=1:size(y,2)
        Psat(i)=exp(a(i)-(b(i)/(T+c(i)))); %kPa
    end       
        
    %solve for phi(i) and x
    T=T+273.15; %K
    B = B_ij(T, Tc, w, Zc, Vc, Pc);
    phi = phi_i(y, B, T, P, Psat);
    for i = 1:size(y,2)
        PHI(i) = (phi(i,1)/phi(i,2));
        x(i) = (y(i)*PHI(i)*P)/(gamma(i)*Psat(i));
    end
    %normalize x
    while sum(x) ~= 1
        for i = 1:size(y,2)
            xnew(i) = x(i)/sum(x);
        end
        x = xnew;
    end
       
    %solve for gamma
    gamma = gamma_UNIFAC(x, T, par, a_mn);        
        
    gammadiff = 1;
    while abs(gammadiff) > tol
        %solve for x
        for i = 1:size(x,2)
            PHI(i) = (phi(i,1)/phi(i,2));
            x(i) = (y(i)*PHI(i)*P)/(gamma(i)*Psat(i));
        end
        while sum(x) ~= 1
            for i = 1:size(y,2)
                xnew(i) = x(i)/sum(x);
            end
            x = xnew;
        end
        gammanew = gamma_UNIFAC(x, T, par, a_mn);
        gammadiff = gammanew-gamma;
        gamma = gammanew;
    end
    
    %calculating Pjsat
    Pjsat=Psat(1);
    for i=1:size(y,2) %1x2
        z(i)=(y(i)*PHI(i)*Pjsat)/(Psat(i)*gamma(i));
    end
    j=sum(z);
    Pjsat=P*j; %kPa
        
    %find the new value of T from Pjsat
    Tnew=((b(1)/(a(1)-log(Pjsat)))-c(1)); %oC
    T=T-273.15;
Tdiff = Tnew-T;
T=Tnew;
disp(T+273.15) %K
int=int+1;
end
disp('No. of iterations =')
disp(int)
d_T = T+273.15;
end
