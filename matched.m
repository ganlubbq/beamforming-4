function [vtM mtD] = matched(mtH,iUsers,iTxAnt)

dBeta = (trace(mtH*mtH'*ones(iUsers,iUsers)))^(-0.5);
vtM = dBeta*mtH'*ones(iUsers,1);

vtM = vtM/norm(vtM);

% rx base
mtD = (1/dBeta)*inv(diag(mtH*mtH'*ones(iUsers,1)));