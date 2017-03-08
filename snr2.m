% dNpow = ((abs(mtD2*mtHro*dDr).^2 + abs(diag(mtD2)).^2)*dNoisePwr);
%  
% snrTwo = (abs(mtD2*mtHr*vtM*dDr).^2) ./ dNpow;

% dNpow = ((abs(mtHro*dDr).^2 + 1)*dNoisePwr);
% 
% snrTwo = (abs(mtHr*vtM*dDr).^2) ./ dNpow;

% snrTwo2 = snrTwo1;
% snrTwo3 = snrTwo1;

function snrTwo = snr2(sys,pa)%mtHro,dDr,vtM,mtHr    sys.mtHro,sys.dDr,sys.M,sys.mtHr

ues = length(sys.mtHro);
snrTwo = zeros(ues,1);

for uu=1:ues
%     dNpow = ((abs(mtD2(uu,uu)*mtHro(uu)*dDr).^2 + abs(mtD2(uu,uu)).^2)*dNoisePwr);
%     snrTwo2(uu) = (abs(mtD2(uu,uu)*mtHr(uu,:)*vtM*dDr).^2) / dNpow;
%     
    dNpow = ((abs(sys.mtHro(uu)*sys.dDr).^2 + 1)*pa.dNoisePwr);
    snrTwo(uu) = (abs(sys.mtHr(uu,:)*sys.M*sys.dDr).^2) / dNpow;
end
% 
% snrTwo - snrTwo1
% snrTwo - snrTwo2
% snrTwo - snrTwo3
% (snrTwo - snrTwo3) ./ snrTwo