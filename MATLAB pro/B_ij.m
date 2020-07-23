function B = B_ij(T, Tc, w, Zc, Vc, Pc)
%B_ij calculates for the values of B which will be used to calculate the
%value of phihat
 
%INPUT:
%T - Temperature in Kelvin
%Tc - row vector of the critical temperature (K) of each species
%w - row vector 
%Zc - row vector
%Vc - row vector of the critical volume (cm3/mol) of each species
%Pc - row vector of the critical pressure (bar) of each species
 
%OUTPUT
%B - square matrix
 
%matrix holders for the values of w_ij, Tc_ij, Zc_ij, Vc_ij, Pc_ij, Tr_ij,
%and B_ij
w_ij = zeros(size(w,2));
Tc_ij = zeros(size(Tc,2));
Zc_ij = zeros(size(Zc,2));
Vc_ij = zeros(size(Vc,2));
Pc_ij = zeros(size(Pc,2));
Tr_ij = zeros(size(Tc,2));
B_ij = zeros(size(Tc,2));
R = 83.14; %unit: (bar cm3/mol K)
%calculates the values of w_ij, Tc_ij, Zc_ij, Vc_ij, Pc_ij, Tr_ij,and B_ij
for i = 1:size(w,2)
    for j = 1:size(w,2)
        w_ij(i,j) = (w(i)+w(j))/2;
        Tc_ij(i,j) = (Tc(i)*Tc(j))^(1/2);
        Zc_ij(i,j) = (Zc(i)+Zc(j))/2;
        Vc_ij(i,j) = (((Vc(i)^(1/3))+(Vc(j)^(1/3)))/2)^3;
        Pc_ij(i,j) = (Zc_ij(i,j)*R*Tc_ij(i,j))/Vc_ij(i,j);
        Tr_ij(i,j) = T/Tc_ij(i,j);
        B0 = 0.083-(0.422/(Tr_ij(i,j)^1.6));
        B1 = 0.139-(0.172/(Tr_ij(i,j)^4.2));
        B_ij(i,j) = ((R*Tc_ij(i,j))/Pc_ij(i,j))*(B0+w_ij(i,j)*B1);
    end
end
B = B_ij;
end
