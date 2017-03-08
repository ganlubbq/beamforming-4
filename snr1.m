% dNpow = (diag(abs(mtD).^2)*dNoisePwr);
% 
% snrOne = (abs(mtD*mtH*vtM).^2) ./ dNpow;

% snrOne1 = (abs(mtH*vtM).^2) ./ dNoisePwr;

function snrOne = snr1(sys,pa)%mtHro,vtM,mtH

ues = length(sys.mtHro);
snrOne = zeros(ues,1);

for uu=1:ues
    snrOne(uu) = (abs(sys.mtH(uu,:)*sys.M).^2) / pa.dNoisePwr;
end