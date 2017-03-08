function chGain = ploss_sha_ant(pDist,shadBool,directBool)
    
angleDbl = (angle(pDist)/pi)*180;        
    
% Antenna gain cf. 3GPP SCM standard having sector 1 as reference
antGainDbl = -min(12*(angleDbl/70)^2, 20);

%chGain = -35.3 - 37.6*log10(abs(pDist)) + 8*randn;
% 3GPP 25814 Dist in KM -128.1
chGain = -35.3 - 37.6*log10(abs(pDist));
if shadBool
    chGain = chGain - 8*randn;
end
if directBool
    chGain = chGain + antGainDbl;
end

chGain = 10 .^ (chGain/10);