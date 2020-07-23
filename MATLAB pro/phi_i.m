function phi = phi_i(y, B, T, P, Psat)
%input P and Psat is in kPa
%input T is in K
%B_ij calculates for the values of B which will be used to calculate the
%value of phihat
 
%INPUT:
%y - row vector of the vapor phase composition
%T - Temperature in Kelvin
%P - Pressure in kPa
%Psat - saturated pressure in kPa
%B - square matrix of the values of B
 
%OUTPUT
%phi - row matrix [phihat phisat]
 
%matrix holders for the values of del_ij, z, phihat, and phisat
del_ij = zeros(size(y,2));
z = zeros(size(y,2));
phihat = zeros(size(y,2),1);
phisat = zeros(size(y,2),1);
R = 8314; %unit: (kPa cm3/mol K)
%calculates the values of del_ij, phihat and phisat
for k = 1:size(y,2)
    for i = 1:size(y,2)
        for j = 1:size(y,2)
            del_ij(i,j) = 2*B(i,j)-B(i,i)-B(j,j);
            z(i,j) = y(i)*y(j)*(2*del_ij(i,k)-del_ij(i,j));
        end
    end
    a = sum(z);
    phihat(k) = exp((P/(R*T))*(B(k,k)+(1/2)*sum(a)));
    phisat(k) = exp((B(k,k)*Psat(k))/(R*T));
end
phi = [phihat phisat];
end
