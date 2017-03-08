function [ vtS ] = MGsinr( pa, sys )
%calculates the minimal SINR for each UE in the MA scheme
    rxPowVtr = zeros(1,pa.iUsers);
    intePowVtr = rxPowVtr;
    for uu=1:pa.iUsers
        gu = sys.ueGrpVtr(uu);
        for gg=1:pa.numOfGroups
            if gg == gu
                rxPowVtr(uu) = (abs(sys.mtH(uu,:)*sys.M{gg}).^2);
            else
                intePowVtr(uu) = intePowVtr(uu) + (abs(sys.mtH(uu,:)*sys.M{gg}).^2);
            end
        end
    end
    snrOne = (rxPowVtr)./((intePowVtr) + pa.dNoisePwr);
    vtS = snrOne;
end

