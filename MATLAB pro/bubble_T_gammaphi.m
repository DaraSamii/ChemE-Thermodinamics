function [b_T,y] = bubble_T_gammaphi(x,Tc,w,Zc,Vc,P,Pc,a,b,c,par,a_mn)
%bubble_T_gammaphi calculates the bubble temperature and vapor phase
%composition of a liquid mixture
 
%INPUT:
%x - row vector of the mole fraction of the components
%P - Pressure in bar
%Tc - row vector of the critical temperature (K) of each species
%w - row vector 
%Zc - row vector
%Vc - row vector of the critical volume (cm3/mol) of each species
%Pc - row vector of the critical pressure (bar) of each species
%a,b,c - row vectors containing the Antoine constants of each species
%par, a_mn - parameters for UNIFAC
 
%OUTPUT
%y - vapor phase composition of the liquid mixture
%b_T - bubble pressure in bar
 
 
%convert the pressure to kPa
P=P*100;
 
%matrix holders for the values of PHI, Psat, z, and y
PHI = ones(size(x,2),1);
Tsat = zeros(size(x,2),1);
Psat = zeros(size(x,2),1);
z = zeros(1,size(x,2));
y = zeros(1,size(x,2));
 
%find the value of Tsat using the given Pressure
for i=1:size(x,2)
    Tsat(i)=((b(i)/(a(i)-log(P)))-c(i));
end
 
%find the initial value of T
for i=1:size(x,2)
    z(i)=((x(i)*Tsat(i)));
end
T=sum(z);
 
%find the value of Psat(i) using Reverse Antoine and
%initial value of T
for i=1:size(x,2)
    Psat(i)=exp(a(i)-(b(i)/(T+c(i))));
end
 
%convert the initial temperature to Kelvin for
%calculation of activity coefficient.
T=T+273.15;
 
%calculate activity coeffucients using UNIFAC and the T initial
gamma = gamma_UNIFAC(x, T, par, a_mn);
 
%calculating Pjsat where j=1
Pjsat=Psat(1);
for i=1:size(x,2)
    z(i)=((x(i)*gamma(i)/PHI(i))*(Psat(i)/Pjsat));
end
j=sum(z);
Pjsat=P/j;
 
%find the new value of T from Pjsat
T=((b(1)/(a(1)-log(Pjsat)))-c(1));
 
%iteration to find the bubble temperature and
%composition of the mixture
 
fprintf('\n\nIteration values of T (bubble T) in Kelvin :\n\n');
tol=10^-6; %tolerance for the iteration
diff=1;
counter=1;
while abs(diff)>tol
    
    %Evaluate Psat and y from the new value of T
    for i=1:size(x,2)
        Psat(i)=exp(a(i)-(b(i)/(T+c(i))));
    end
    
    %solve for the vapor phase composition of the mixture
    for i=1:size(x,2)
        y(i)=((x(i)*gamma(i)*Psat(i))/(PHI(i)*P));
    end
    
    T=T+273.15; %convert the temperate to Kelvin
    
    %calculate the value of B and phi(hat) and phi(sat)
    %from subroutine functions
    B = B_ij(T, Tc, w, Zc, Vc, Pc);
    phi = phi_i(y, B, T, P, Psat);
    
    %calculate activity coefficients from UNIFAC
    for i = 1:size(x,2)
        PHI(i) = (phi(i,1)/phi(i,2));
    end
    gamma = gamma_UNIFAC(x, T, par, a_mn);
    
    %calculating Pjsat from the previous activity coefficient
    Pjsat=Psat(1);
    for i=1:size(x,2)
        z(i)=((x(i)*gamma(i)/PHI(i))*(Psat(i)/Pjsat));
    end
    j=sum(z);
    Pjsat=P/j;
 
    %find the new value of T from Pjsat
    Tnew=((b(1)/(a(1)-log(Pjsat)))-c(1));
    T=T-273.15;
    disp(T+273.15);
    counter=counter+1;
    %check for the tolerance criterion
    diff=(T-Tnew);
    T=Tnew;
    
end
disp('No of iterations: ');
disp(counter);
b_T = T+273.15;
end
