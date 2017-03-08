function [vtM mtD] = zforcing(mtH,iUsers,iTxAnt)
% Zero forcing
dBeta = (trace(inv(mtH*mtH')*ones(iUsers,iUsers)))^(-0.5);
vtM = dBeta*mtH'*inv(mtH*mtH')*ones(iUsers,1);
mtD = (1/dBeta).*eye(iUsers);