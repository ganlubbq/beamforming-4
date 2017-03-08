function [ new_filter ] = limitRelayFilterMulti( pa, sys )
    % depends ont the power of sys.M
    totalPower = 0;
    for gg=1:pa.numOfGroups
        totalPower = totalPower + norm(sys.MD*sys.Hbr*sys.M{gg})^2;
    end        

    totalPower = totalPower + pa.dNoisePwr*norm(sys.MD)^2;
    
    new_filter = sqrt(pa.dRelayPower/totalPower)*sys.MD;    
end

