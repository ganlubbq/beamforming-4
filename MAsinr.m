function [ vtS ] = MAsinr( pa, sys )
%calculates the minimal SINR for each UE in the MA scheme
%     rxPowVtr = zeros(1,pa.iUsers);
%     intePowVtr = rxPowVtr;
%     rxPow2Vtr = rxPowVtr;
%     intePow2Vtr = rxPowVtr;
%     for uu=1:pa.iUsers
%         gu = sys.ueGrpVtr(uu);
%         for gg=1:pa.numOfGroups
%             if gg == gu
%                 rxPowVtr(uu) = (abs(sys.mtH(uu,:)*sys.M{gg}).^2);
%                 rxPow2Vtr(uu) = (abs(sys.Hru(uu,:)*sys.MD*sys.Hbr*sys.M{gg}).^2);
%             else
%                 intePowVtr(uu) = intePowVtr(uu) + (abs(sys.mtH(uu,:)*sys.M{gg}).^2);
%                 intePow2Vtr(uu) = intePow2Vtr(uu) + (abs(sys.Hru(uu,:)*sys.MD*sys.Hbr*sys.M{gg}).^2);
%             end
%         end
%         intePow2Vtr(uu) = intePow2Vtr(uu) + norm(sys.Hru(uu,:)*sys.MD)^2*pa.dNoisePwr;
%     end
%     snrOne = (rxPowVtr*pa.dBasePower)./((intePowVtr*pa.dBasePower) + pa.dNoisePwr);
%     snrTwo = (rxPow2Vtr*pa.dBasePower)./((intePow2Vtr*pa.dBasePower) + pa.dNoisePwr);
%     vtS = snrOne + snrTwo;
%% TODO: A POTENCIA JA TEM QUE ESTAR NO FILTRO
    rxPowVtr = zeros(1,pa.iUsers);
    intePowVtr = rxPowVtr;
    rxPow2Vtr = rxPowVtr;
    intePow2Vtr = rxPowVtr;
    for uu=1:pa.iUsers
        gu = sys.ueGrpVtr(uu);
        for gg=1:pa.numOfGroups
            if gg == gu
                rxPowVtr(uu) = (abs(sys.mtH(uu,:)*sys.M{gg}).^2);
                rxPow2Vtr(uu) = (abs(sys.Hru(uu,:)*sys.MD*sys.Hbr*sys.M{gg}).^2);
            else
                intePowVtr(uu) = intePowVtr(uu) + (abs(sys.mtH(uu,:)*sys.M{gg}).^2);
                intePow2Vtr(uu) = intePow2Vtr(uu) + (abs(sys.Hru(uu,:)*sys.MD*sys.Hbr*sys.M{gg}).^2);
            end
        end
        % add the transmit powers powers
        %rxPowVtr(uu) = rxPowVtr(uu)*pa.dBasePower;
        %rxPow2Vtr(uu) = rxPow2Vtr(uu)*pa.dBasePower;
        %intePowVtr(uu) = intePow2Vtr(uu)*pa.dBasePower;
        %intePow2Vtr(uu) = intePow2Vtr(uu)*pa.dBasePower;
        %
        intePow2Vtr(uu) = intePow2Vtr(uu) + norm(sys.Hru(uu,:)*sys.MD)^2*pa.dNoisePwr;
    end
    snrOne = (rxPowVtr)./((intePowVtr) + pa.dNoisePwr);
    snrTwo = (rxPow2Vtr)./((intePow2Vtr) + pa.dNoisePwr);
    vtS = snrOne + snrTwo;
end