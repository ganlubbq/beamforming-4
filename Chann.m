classdef Chann
    properties
        mtH
        mtHro
        mtHir
        M
        D
        D2
        mtHr
        Hbr
        Hru
        dDr
        MD
        stPos
        ueGrpVtr
    end
    methods
        function sys = Chann(pa)
            %% Generates channel
            
            %pa = cheat.pa();
            % defines the positions of the base, relay and users
            [sys.stPos] = placesNodes(pa.dR,pa.dHotRadius,pa.dR/2,pa.iUsers,pa.hotUE,pa.posScheme);
            
            % calculate path loss and shadowing
            vtPlossBase = zeros(1,pa.iUsers);
            vtPlossRelay = vtPlossBase;
            for uu=1:pa.iUsers
                vtPlossBase(uu) = ploss_sha_ant(sys.stPos.vtPosUsers(uu) - sys.stPos.cpPosBase,true,true);
                vtPlossRelay(uu) = ploss_sha_ant(sys.stPos.vtPosUsers(uu) - sys.stPos.cpPosRelay,true,false);
            end
            
            dPlossRelay = ploss_sha_ant(sys.stPos.cpPosRelay - sys.stPos.cpPosBase,true,true);
            
            % Fading channel realization
            sys.mtH = los_channel(pa.iUsers, pa.iTxAnt, 1e-5, 2*pi/3); %los_channel(iUsers, iTxAnt, dK, 2*pi/3);
            
            %TODO: change single relay antenna algorithms to complain
            %abount relay channels with wrong sizes
            %sys.mtHir = los_channel(pa.relayTxAnt, pa.iTxAnt, 1e-5, 2*pi/3);            
            %sys.mtHro = los_channel(pa.iUsers, pa.relayTxAnt, 1e-5, 2*pi/3);
            sys.Hbr = los_channel(pa.relayTxAnt, pa.iTxAnt, 1e-5, 2*pi/3);
            sys.Hru = los_channel(pa.iUsers, pa.relayTxAnt, 1e-5, 2*pi/3);
            
            % aplies the PL and shadowing to the channels
            sys.Hbr = sys.Hbr*dPlossRelay;
            
            for uu=1:pa.iUsers
                sys.Hru(uu,:) = sys.Hru(uu,:)*vtPlossRelay(uu);
                sys.mtH(uu,:) = sys.mtH(uu,:)*vtPlossBase(uu);
            end
            
            sys.mtHir = sys.Hbr(1,:);            
            sys.mtHro = sys.Hru(:,1);
            
            % tx rx relay
            % TODO: remove the general matrix and put it in the specific
            % algorithms. This does not make sense with a multiple antenna
            % relay
            sys.mtHr = sys.mtHro * sys.mtHir; % h_3 * h_2
            
            %sys.lH = sys.mtH';
            
            %% fixed relay gain
            % the gaim dependes in the chosen transmit filter according to
            % sqrt(pa.dRelayPower/(pa.dBasePower*norm(sys.mtHir*vtM)^2 +
            % pa.dNoisePwr));
            %% TODO: calculade the correct dr
            sys.dDr = NaN; %sqrt(pa.dRelayPower/(pa.dBasePower*norm(sys.mtHir)^2 + pa.dNoisePwr));
            sys.MD = NaN;
            disp('Channel generated.');
        end
    end
end

function mtH = los_channel(iN, iM, dK, dAngleRange) % , iPL, dR

dDist = 0.5;
vtElements = 0:(iM-1);
vtAngles = dAngleRange*rand(iN,1);
mtH = zeros(iN, iM);

for i=1:iN,
    mtH(i,:) = exp(sqrt(-1)*2*pi*dDist*cos(vtAngles(i))*vtElements);
end

mtHw = sqrt(0.5).*(randn(iN, iM) + randn(iN, iM).*sqrt(-1));

if dK > 1e6
    mtHw = zeros(iN,iM);
end

mtH = sqrt(dK/(1+dK))*mtH + sqrt(1/(1+dK))*mtHw;
end