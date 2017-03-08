%% stores the simulation parameters to be used throught the simulation
classdef Parameters < handle
    properties
        iIter;
        dR;
        iTxAnt;
        simulation;
        iUsers;
        subsim;
        iSymbols;
        getBER;
        minSINRdB;
        posScheme;
        hotUEperc;
        hotUE;
        dHotRadius;
        numOfGroups;
        vtScheme;
        vtAlgs;
        % 0, Matched filter
        % 2, Zero forcing
        % 5, Switched fixed beams
        % 6, USMF
        % 8, SDR method
        % 9, Rand
        % 10, Iterative update with RC
        % 11, Matched filter for worst user
        dNoisePwr;
        dNoiseStd;
        % linksim
        dBasePower;
        dRelayPower;
        relayTxAnt;
    end
    methods
        function pa = Parameters()
            pa.iIter = 3; %1e2; %2.5e2;
            pa.dR = 0.334;
            pa.iTxAnt = 3;
            pa.relayTxAnt = 2;
            pa.simulation = 'test';
            pa.iUsers = 6;
            pa.subsim = '001';
            
            %% Definition of fixed parameters
            %pa.iN = 4;
            %pa.sModulation = 'PSK';
            pa.iSymbols = 1; %100;                        
            %pa.sKstr = 'NLOS'; %{'NLOS','LOS'};
            %pa.dK = 1e-5; % [ 1e-5 1e5 ];
            % ber results take aprox. 10x more time
            pa.getBER = false;
            % fFixedSeed = false;
            pa.minSINRdB = 10;
            pa.posScheme = 'random_hex'; %{'circleOnRelay','random_circ','random_hex'};
            
            %% hotspot only
            pa.hotUEperc = 0.45;
            pa.hotUE = round(pa.hotUEperc*pa.iUsers);
            pa.dHotRadius = 0.3*pa.dR;
            
            %% multigroup only
            pa.numOfGroups = 3;
            
            if not(pa.getBER)
                pa.iSymbols = 1;
            end
            
            % Noise for one LTE RB (180kHz)
            % noise figure 9 dB
            pa.dNoisePwr = 10^((-174 + 9)/10.0 - 3)*180e3;
            pa.dNoiseStd = sqrt(pa.dNoisePwr);
            
            pa.vtScheme = {'MGsdr','MGBh','MGRIR'};%{'rotatesdr','noRelay','rotate','iterative','MRC'};
            %{'novo','noRelay','MRC','iterative'}; %{'novo','noRelay','onlyRelay','eqCombine','MRC','iterative'};
            pa.vtAlgs = [8 10];%6 [0 6 8 10];
            disp('Parameters generated.');
        end
    end
end


