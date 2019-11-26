function [ cij, depp ] = ACFLNGBHE_to_cij( A,C,F,L,N,Gc,Bc,Hc,Ec,depths,dep )
% JBR 12/07/17
% Calculate cij from velocity parameters
%
cij = zeros(6,6);
[val, idep] = min(abs(depths-dep));
depp = depths(idep);

cij(1,1) = A(idep) + Bc(idep) + Ec(idep);
cij(2,2) = A(idep) - Bc(idep) + Ec(idep);
cij(3,3) = C(idep);
cij(4,4) = L(idep) - Gc(idep);
cij(5,5) = L(idep) + Gc(idep);
cij(6,6) = N(idep) - Ec(idep);
cij(1,2) = A(idep) - 2*N(idep) - Ec(idep);
cij(1,3) = F(idep) + Hc(idep);
cij(2,3) = F(idep) - Hc(idep);

cij(2,1) = cij(1,2);
cij(3,1) = cij(1,3);
cij(3,2) = cij(2,3);

end

