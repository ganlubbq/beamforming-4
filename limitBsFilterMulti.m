function [ new_filter ] = limitBsFilterMulti( pa, sys )
    current_filter = sys.M;

    new_filter = current_filter;
    totalPower = 0;
    % generate random tx vectors
    for gg=1:pa.numOfGroups
        totalPower = totalPower + norm(current_filter{gg})^2;
    end

    % limit tx power
    for gg=1:pa.numOfGroups
        new_filter{gg} = current_filter{gg}*sqrt(pa.dBasePower/totalPower);
    end
                
end

