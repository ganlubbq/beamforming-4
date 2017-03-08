function [vtM, mtD] = doPrecoding(iAlg,mtH,iUsers,iTxAnt)

sAlgName= algName(iAlg);

[vtM mtD] = feval(sAlgName,mtH,iUsers,iTxAnt);
