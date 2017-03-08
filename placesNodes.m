function [stPos] = placesNodes(dR,dHotRadius,dRelayDist,iNumUser,iNumHotUser,sMode)

% This funcion works in a sector intead of a cell.
% Thus, it places users inside of an third of exagon with the base station in its
% left border.

% dR is the exagon diameter

stPos.cpPosBase = complex(0,0);
stPos.cpPosRelay = complex(dRelayDist,0);
if strcmp(sMode,'circleOnRelay')
    dAng = (2*pi)/iNumUser;
    vtAux = 0:(iNumUser-1);
    vtAux = exp(dAng*1i.*vtAux);
    dUsrR = dR - dRelayDist;
    stPos.vtPosUsers = dUsrR*vtAux + stPos.cpPosRelay;
elseif strcmp(sMode,'random_circ')
    for uu=1:iNumUser
        retry = true;
        while retry
            pX = rand*dR - dR/2;
            pY = rand*dR - dR/2;
            radi = sqrt(pX^2 + pY^2);
            if radi <= dR/2
                retry = false;
            end
        end
        stPos.vtPosUsers(uu) = pX + dR/2 + 1i*pY;
    end
elseif strcmp(sMode,'random_hex')
    for uu=1:iNumUser
        validPosBool = false;
        % Sort position and test if it's within the hexagon
        while validPosBool == false,
            posCpx = 2*dR*complex(rand-0.5, rand-0.5);
            
            angleDbl = angle(posCpx);
            distanceDbl = abs(posCpx);
            limitDbl = hexBorder(dR, angleDbl);
            
            if(distanceDbl < limitDbl)
                validPosBool = true;
            end
        end
        stPos.vtPosUsers(uu) = posCpx;
    end
    
    stPos.vtPosUsers = ((stPos.vtPosUsers + dR)/2);
    
elseif strcmp(sMode,'random_hex_w_hot')
    genUsers = iNumUser - iNumHotUser;
    
    % pesitions in general ues a hexagon
    for uu=1:genUsers
        validPosBool = false;
        % Sort position and test if it's within the hexagon
        while validPosBool == false,
            posCpx = 2*dR*complex(rand-0.5, rand-0.5);
            
            angleDbl = angle(posCpx);
            distanceDbl = abs(posCpx);
            limitDbl = hexBorder(dR, angleDbl);
            
            if(distanceDbl < limitDbl)
                validPosBool = true;
            end
        end
        stPos.vtPosUsers(uu) = posCpx;
    end
    
    % pesitions in hotspot ues a hexagon
    for uu=1:iNumHotUser
        validPosBool = false;
        % Sort position and test if it's within the hexagon
        while validPosBool == false,
            posCpx = 2*dHotRadius*complex(rand-0.5, rand-0.5);
            
            angleDbl = angle(posCpx);
            distanceDbl = abs(posCpx);
            limitDbl = hexBorder(dHotRadius, angleDbl);
            
            if(distanceDbl < limitDbl)
                validPosBool = true;
            end
        end
        stPos.vtPosUsers(genUsers + uu) = posCpx;
    end
    
    % center general ues
    stPos.vtPosUsers(1:genUsers) = ((stPos.vtPosUsers(1:genUsers) + dR)/2);
    
    % center hotspot ues
    stPos.vtPosUsers((genUsers+1):end) = (stPos.vtPosUsers((genUsers+1):end)/2 + dRelayDist);
    
else
    disp('Invalid placement mode!');
    exit
end

% hold on
% plot(stPos.cpPosBase,'o');
% plot(stPos.cpPosRelay,'s');
% plot(stPos.vtPosUsers,'x');
% hold off
% axis equal
end

function vtR = hexBorder(dRadius, vtTheta)

for i=1:length(vtTheta),
    dTheta = vtTheta(i);
    if dTheta < 0,
        dTheta = 2*pi + dTheta;
    end
    if dTheta > pi/3,
        dTheta = mod(dTheta, pi/3);
    end
    
    vtR(i) = dRadius / (cos(dTheta) + sin(dTheta)/tan(pi/3));
end
end