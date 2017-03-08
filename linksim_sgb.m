%% Function linksim_sgb
% Link level simulator for single-group multicast beamforming
% Input parameters:
%   iN:             Modulation order
%   sModulation:    String containing the modulation name ('PSK' or 'QAM')
%                   Note: the tau parameter of THP is only implemented for
%                   QPSK and 16-QAM, for other modulation schemes vtTau
%                   needs to be properly adjusted
%   dSNR:           Eb/No in dB
%   iSymbols:       Number of symbols during which the channel can be assumed invariant
%   iIter:          Number of channel realizations to be simulated
%   iTxAnt:         Number of transmit antennas at the base station
%   iUsers:         Total number of single-antenna users
%   iAlg:           Identifies the selected single-group beamforming algorithm
%                   0- Matched filter, 1- Matched average filter, 2- Zero-forcing,
%                   3- MMSE, 4- ZF-THP, 5- Switched fixed beams, 6- USMF,
%                   7- USMF with correlation, 8- SDR method, 9- Pure randomization
%   dK:             Ricean fading parameter
%   iPL:            0- Without path-loss, 1- With path-loss
%   dR:             Cell radius
%   bBER            BER results flag
%   iRelayScheme    String containing the relay scheme (noRelay,onlyRelay,eqCombine,MRC)

function [stRes, cellAux] = linksim_sgb(pa, iAlg, iRelayScheme)

%% read parameters
cellAux = cell(1,0);
%iN = pa.iN;
%sModulation = pa.sModulation;
iSymbols = pa.iSymbols;
%iIter = pa.iIter;
iTxAnt = pa.iTxAnt;
iUsers = pa.iUsers;
%dK = pa.dK;
bBER = pa.getBER;

%stChann =  struct(mtH,[],mtHro,[],mtHir,[],vtM,[],mtD,[],mtD2,[],mtHr,[]);

% Calculate average energy per symbol for the default constellation
% 36.104 3GPP max 18 dbm (4 ant.)
% 25814 43 dbm Tot 5 MHz / 6 RBs
%dBasePower = 19.9526/25;
%dBasePower = (10^((43 - 30)/10))/12;

%% Calculates tx power for minimum SINR ate the border (considering Ploss and sectorization)
%% TODO: works with random_hex

pa.dBasePower = (pa.minSINRdB*pa.dNoisePwr)/ploss_sha_ant(pa.dR,false,true);

%% TODO: relay power must obey an ratio with base power and it must change
% in a loop
%dRelayPower = dBasePower;

% currently designed for the same rule as the base (using the hotspot
% radius)
pa.dRelayPower = (pa.minSINRdB*pa.dNoisePwr)/ploss_sha_ant(pa.dHotRadius,false,true);

% puts the power back into the base for fair comparisson
if strcmp(iRelayScheme,'SnoRelay')
    pa.dBasePower = pa.dBasePower + pa.dRelayPower;
    pa.dRelayPower = 0;
end

% 500 m ISD -> radium 334 m

%% calculate sysmbol energy
% Generate constellation
% if strcmp(sModulation, 'QAM') == 1,
%     vtConst = qammod(0:iN-1, iN, 0, 'gray').';
%     dEsAvg = sum(abs(vtConst).^2) / iN;
%     vtConst = vtConst./sqrt(dEsAvg);
% elseif strcmp(sModulation, 'PSK') == 1,
%     vtConst = pskmod(0:iN-1, iN, pi/4, 'gray').';
% end
% vtConst = vtConst * sqrt(pa.dBasePower);
% dEsAvg = sum(abs(vtConst).^2) / iN;
% pa.dEsAvg = dEsAvg;

%% Noise for one LTE RB (180kHz)
% noise figure 9 dB
%dNoisePwr = pa.dNoisePwr;
%dNoiseStd = pa.dNoiseStd;

% Subcarrier noise power per subcarrier (in Watts)
% Noise power density in dBm + Noise Figure in dB converted to linear scale
% times the bandwidth of one subcarrier
% TODO: WHERE DO THESE VALUES COME FROM?
% paramSt.noisePwrDbl = 10^((-174 + 9)/10.0 - 3)*paramSt.PRBBandInt;

%% Initialize output variables
if bBER
    vtBER = zeros(pa.iIter,1);
    vtMaxBER = zeros(pa.iIter,1);
end
mtSNR = zeros(pa.iIter,iSymbols,iUsers);

if strcmp(iRelayScheme,'iterative') == 1,
    cellAux{1}{1} = zeros(1,pa.iIter);
end

if ~isempty(regexp(iRelayScheme, 'rotate', 'once'))
    cellAux{1}{1} = zeros(1,pa.iIter);
end

%% generates the multiple random streams. One for each iteration
mainStream = RandStream.getDefaultStream;
streamVector = cell(1,pa.iIter);
for i=1:pa.iIter,
    streamVector{i} = RandStream('mt19937ar','seed',sum(100*rand));
    %streamVector{i}.rand(1,3)
end

%% channel realization loop
for i=1:pa.iIter,
    
    %% display iteration information
    %if mod(i,1000) == 0
    disp(['Iter: ' num2str(i)]);
    %end
    
    %% set the random number generator state for this iteration
    RandStream.setDefaultStream(streamVector{i});
    
    %xxxxxxxxxxxxxxxxxxxx
    sys = Chann(pa);  %cheat.sys();
    
    %     if bBER
    %         % Get random constellation symbols
    %         vtXsi = randint(1,iSymbols,iN);
    %         vtXsi = repmat(vtXsi, iUsers, 1);
    %         vtX = vtConst(vtXsi+1);
    %
    %
    %         % Generate noise
    %         vtNoise1 = dNoiseStd.*sqrt(0.5).*(randn(iUsers,iSymbols) + randn(iUsers,iSymbols).*sqrt(-1));
    %         vtNoise2 = dNoiseStd.*sqrt(0.5).*(randn(iUsers,iSymbols) + randn(iUsers,iSymbols).*sqrt(-1));
    %         dNoiseR = dNoiseStd.*sqrt(0.5).*(randn(iUsers,iSymbols) + randn(iUsers,iSymbols).*sqrt(-1));
    %         %if bRelay
    %         %    vtRNoise = dRNoiseStd.*sqrt(0.5).*(randn(1,iSymbols) + randn(1,iSymbols).*sqrt(-1));
    %         %end
    %
    %         % Reset error count
    %
    %         iErrors = 0;
    %         vtErrors = zeros(iUsers,1);
    %     end
    
    %% TODO: remove power normalization
    
    %% Tx/Rx filter generation
    
    % currently the filter generation is divided in two case:
    % -the general MA case were no filter normalization is used
    % this is the correct and more general sice the SINR cant be calculated
    % using an normalized filter in all case
    % -legacy normalizated. All previously coded solutions that already
    % used normalized filters
    if (~isempty(strmatch('MA',iRelayScheme)))
        % multi antenna relayed multigroup mode
        
        %% seperar grupos
        % calcula UEs por grupo
        grpUes = floor(pa.iUsers/pa.numOfGroups);
        uesRemaining = pa.iUsers - grpUes*pa.numOfGroups;
        
        %vetor com o numero do grupo de cada UE
        sys.ueGrpVtr = zeros(1,pa.iUsers);
        for gg=1:pa.numOfGroups
            sys.ueGrpVtr((1:grpUes) + (gg-1)*grpUes) = gg;
        end
        sys.ueGrpVtr((end-uesRemaining):end) = pa.numOfGroups;
        
        %% algoritmos
        if strcmp(iRelayScheme,'MArand') == 1,
            
            %[sys.M, sys.MD] = MArand(sys,pa);
            sys.M = cell(1,pa.numOfGroups);
            
            % generate random tx vectors
            for gg=1:pa.numOfGroups
                sys.M{gg} = rand(pa.iTxAnt,1) + 1i*rand(pa.iTxAnt,1);
            end
                                    
            % generate random relay matrix
            sys.MD = rand(pa.relayTxAnt) + 1i*rand(pa.relayTxAnt);                       
            
            %elseif strcmp(iRelayScheme,'MAbrute') == 1,
            %[sys.M, sys.MD] = Brute(sys,pa);
        elseif strcmp(iRelayScheme,'MAga') == 1,
            [sys.M, sys.MD] = MAga(sys,pa);
        end
        % limit tx power
        [ sys.M ] = limitBsFilterMulti( pa, sys );
        
        % limit relay power
        [ sys.MD ] = limitRelayFilterMulti( pa, sys );
    else
        if strcmp(iRelayScheme,'SnoRelay') == 1, % without relay
            
            sys.mtHr = sys.mtHro * sys.mtHir;
            sys.D2 = zeros(iUsers);
            sys.dDr = 0;
            
            [sys.M, sys.D] = doPrecoding(iAlg,sys.mtH,iUsers,iTxAnt);
            
        else % with realy
            sAlgName= algName(iAlg);
            
            if strcmp(iRelayScheme,'SonlyRelay') == 1,
                
                sys.D = zeros(iUsers);
                
                sys.mtHr = sys.mtHro * sys.mtHir;
                
                [sys.M, sys.D2] = feval(sAlgName,sys.mtHr,iUsers,iTxAnt);
                
                sys.mtHef = mtD*sys.mtH + sys.D2*sys.mtHr*sys.dDr;
                
            elseif strcmp(iRelayScheme,'SMRC') == 1,
                [sys.M, sys.D] = feval(sAlgName,sys.mtH,iUsers,iTxAnt);
                
                % checks if relay power is less than or equal to dRelayPower
                pr = pa.dBasePower*abs(sys.dDr*sys.mtHir*sys.M).^2 + sys.dDr*sys.dDr*pa.dNoisePwr;
                if pr > dRelayPower
                    keyboard;
                end
                
                sys.mtHr = sys.mtHro * sys.mtHir;
                
                [sys.D, sys.D2] = mrc(sys.mtHr,sys.mtHro,sys.mtH,sys.dDr,sys.M,pa.dNoisePwr);
                
                sys.mtHef = sys.D*sys.mtH + sys.D2*sys.mtHr*sys.dDr;
                
            elseif strcmp(iRelayScheme,'Siterative') == 1,
                
                % new iterative
                iterSche;
                
            elseif strcmp(iRelayScheme,'rotate') == 1,
                
                sys.dDr = sqrt(pa.dRelayPower/(pa.dBasePower*norm(sys.mtHir)^2 + pa.dNoisePwr));
                % complex rotation
                rotSche;
                
            elseif strcmp(iRelayScheme,'rotate2') == 1,
                
                % complex rotation
                rotSche2;
                
            elseif strcmp(iRelayScheme,'rotatesdr') == 1,
                
                % ALGORITMO ASSUME UM GANHO DE RELAY SUBOTIMO (DE PIOR CASO)
                % QUE NÃO DEPENDIA DO BEAMFORMER
                keyboard;
                
                h = cell(1,iUsers);
                %hx = h;
                g = h;
                %gg = g;
                a = abs(max(max(sys.mtH))) + abs(max(max(sys.mtHr)));
                for uu=1:iUsers
                    %h{uu} = dEsAvg*sys.mtH(uu,:).' / sqrt(pa.dNoisePwr);
                    %hx{uu} = dEsAvg*sys.mtHr(uu,:).' *dDr / sqrt((abs(sys.mtHro(uu)*dDr).^2 + 1)*pa.dNoisePwr);
                    %g{uu} = h{uu}*h{uu}' + hx{uu}*hx{uu}';
                    g{uu} = (1/a)*((pa.dBasePower/pa.dNoisePwr)*(conj(sys.mtH(uu,:).')*(sys.mtH(uu,:)) + ((sys.dDr^2)/(abs(sys.mtHro(uu)*sys.dDr).^2 + 1))*(conj(sys.mtHr(uu,:).')*(sys.mtHr(uu,:)))));
                    %                 disp(max(eig(g{uu})));
                end
                
                %             if i == 364
                %                 keyboard;
                %             end
                
                [sys.M, bRandom] = sdrrel(g,iUsers,iTxAnt);
                
                [sys.D, sys.D2] = mrc(sys.mtHr,sys.mtHro,sys.mtH,sys.dDr,sys.M,pa.dNoisePwr);
                
                %             vtSU = zeros(1,iUsers);
                %
                %             for uu=1:iUsers
                %                 vtSU(uu) = trace((vtM*vtM')*g{uu});
                %             end
                
                cellAux{1}{1}(i) = bRandom;
                
            elseif (~isempty(strmatch('MG',iRelayScheme)))
                %% Multigrupo
                
                % simple multigroup does not use relay
                sys.dDr = 0;
                
                %% seperar grupos
                % calcula UEs por grupo
                grpUes = floor(pa.iUsers/pa.numOfGroups);
                uesRemaining = pa.iUsers - grpUes*pa.numOfGroups;
                
                %vetor com o numero do grupo de cada UE
                sys.ueGrpVtr = zeros(1,pa.iUsers);
                for gg=1:pa.numOfGroups
                    sys.ueGrpVtr((1:grpUes) + (gg-1)*grpUes) = gg;
                end
                sys.ueGrpVtr((end-uesRemaining):end) = pa.numOfGroups;
                
                if strcmp(iRelayScheme,'MGBh') == 1,
                    [sys.M] = Bornhorst(sys,pa);
                    
                elseif strcmp(iRelayScheme,'MGsdr') == 1,
                    %% alg MGsdr
                    
                    [sys.M] = sdrmulti(sys,pa);                                        
                    
                    %% TODO: cellD{gg} = inv(diag(mtH*cellM{gg}));
                    
                    %                 vtM = cell(1,pa.numOfGroups);
                    %                 mtD = vtM;
                    %                 mtD2 = vtM;
                    %                 for gg=1:pa.numOfGroups
                    %                     vtM{gg} = rand(pa.iTxAnt,1) + 1i*rand(pa.iTxAnt,1);
                    %                     vtM{gg} = vtM{gg}/norm(vtM{gg});
                    %
                    %                     selUesVtr = (ueGrpVtr == gg);
                    %
                    %                     [mtD{gg} mtD2{gg}] = mrc(sys.mtHr(selUesVtr,:),sys.mtHro(selUesVtr,:),sys.mtH(selUesVtr,:),dDr,vtM{gg},pa.dNoisePwr);
                    %                 end
                    
                elseif strcmp(iRelayScheme,'MGRIR') == 1,
                    %% ideia de solução
                    % NÃO FUNCIONA DIREITO.
                    % DESEMPENHO RUIM
                    
                    % 3 passos ( não necessriamanete todos caontecem o mesmo
                    % numero de vezers)
                    % 1-maximizar o sinal util alinhando o filtro com os piores
                    % UEs
                    % 2-remover a interferencia de um dado grupo removendo
                    % parte da parte que não projeta no espaçõ nulo dos piores
                    % UEs do filtro do outro grupo
                    % 3-Ajusta de faze da constante compleza que multiplica
                    % cada filtro. Não altera o sinal util mas tedus a
                    % interferencia entre grupos
                    
                    % alg MGRIR
                    [vtM] = mgrir(sys,pa);
                elseif strcmp(iRelayScheme,'MGga') == 1,
                    [sys.M] = MGga(sys,pa);
                end
                
                %cellAux{1}{1}(i) = bRandom;
                % limit tx power
                [ sys.M ] = limitBsFilterMulti( pa, sys );
                    
                % limit relay power
                %[ sys.MD ] = limitRelayFilterMulti( pa, sys );
                
            elseif strcmp(iRelayScheme,'test') == 1,
                rand(1,10)
%                 sys.mtHr = sys.mtHro * sys.mtHir;
%                 
                 sys.M = zeros(iTxAnt,1);
                 %mtD = eye(iUsers);
%                 [vtM1, mtD1] = doPrecoding(8,sys.mtH,iUsers,iTxAnt);
%                 vtSNRi1 = (pa.dBasePower*abs(mtD1*sys.mtH*vtM1).^2) ./ (pa.dNoisePwr*abs(diag(mtD1)).^2);
%                 m1 = min(vtSNRi1);
%                 [vtM2, mtD2] = doPrecoding(10,sys.mtH,iUsers,iTxAnt);
%                 vtSNRi2 = (pa.dBasePower*abs(mtD2*sys.mtH*vtM2).^2) ./ (pa.dNoisePwr*abs(diag(mtD2)).^2);
%                 m2 = min(vtSNRi2);
%                 
%                 if m2 < m1
%                     keyboard;
%                 end
%                 
%                 vtM = vtM1;
%                 mtD = mtD1;
%                 %[vtM, mtD] = doPrecoding(iAlg,sys.mtH,iUsers,iTxAnt);
%                 
%                 
                 %mtD2 = zeros(iUsers);
                 sys.dDr = 0;
%                 %disp(num2str(rand))
            else
                disp([iRelayScheme ' scheme not implemented!'])
                exit
            end
        end
        
        % garantees that the filter is normalized
%         if (iscell(sys.M))
%             [ sys.M ] = normalizeBsFilterMulti( pa, sys );
%         else
%             sys.M = sys.M/norm(sys.M);
%         end
    end
    
    %vtSNRi = zeros(iUsers,1);
    
    %% Transmit each symbol
    for j=1:iSymbols,
        % Tx/Rx processing
        
        %         if strcmp(iRelayScheme,'noRelay') == 1, % without relay
        %             if iAlg == 4,
        %                 % THP
        %                 vtV = zeros(iUsers,1);
        %                 vtXo = mtOrder*vtX(:,j);
        %                 for k=1:iUsers,
        %                     mtProd = mtF*vtV;
        %                     vtV(k) = vtXo(k) + mtProd(k);
        %                     vtV(k) = thp_mod(vtV(k), dTau);
        %                 end
        %                 dBeta = sqrt(pa.dBasePower / trace(diag(diag(mtL).^(-2))*vtV*vtV'));
        %                 mtM = dBeta*mtMn;
        %                 mtD = (1/dBeta).*mtOrder*eye(iUsers);
        %
        %                 vtY = mtD*sys.mtH*mtM*vtV + mtD*vtNoise1(:,j);
        %                 vtY = thp_mod(vtY, dTau);
        %                 vtSNRi = abs(mtD*sys.mtH*mtM*vtV).^2 ./ abs(mtD*vtNoise(:,j)).^2;
        %             else
        %                 if bBER
        %                     vtY = mtD*sys.mtH*vtM*vtX(1,j) + mtD*vtNoise1(:,j);
        %                 end
        %                 vtSNRi = (pa.dBasePower*abs(mtD*sys.mtH*vtM).^2) ./ (pa.dNoisePwr*abs(diag(mtD)).^2);
        %             end
        %         else % with realy
        
        %         if bBER
        %             % TODO: modificar para o novo modelo de ruidos
        %             vtY = sys.mtHef*vtM*vtX(1,j) + vtNoise1(:,j);
        %         end
        %         %vtSNRi = abs(sys.mtHef*vtM*vtX(1,j)).^2 ./ abs(vtNoise1(:,j)).^2;
        %         dNpowT = ((abs(mtD2*sys.mtHro*dDr).^2 + abs(diag(mtD)).^2 + abs(diag(mtD2)).^2)*pa.dNoisePwr);
        %         vtSNRi = (abs(mtD*sys.mtH*vtM + mtD2*sys.mtHr*vtM*dDr).^2)*pa.dBasePower;
        %         vtSNRi = vtSNRi ./ dNpowT;
        
        %% test
        %         acor = zeros(iUsers,pa.numOfGroups);
        %         npcor = acor;
        %         for uu=1:iUsers
        %             %gu = sys.ueGrpVtr(uu);
        %             %selUesVtr = (sys.ueGrpVtr == gg);
        %             for gg=1:pa.numOfGroups
        %                 acor(uu,gg) = abs(sys.mtH(uu,:)*vtM{gg})^2;
        %                 npcor(uu,gg) = abs(sys.mtH(uu,:)*vtM{gg}/(norm(sys.mtH(uu,:))*norm(vtM{gg})))^2;
        %             end
        %         end
        %         acor
        %         npcor
        %
        %         acor = zeros(iUsers,iUsers);
        %         for uu=1:iUsers
        %             for ua=1:iUsers
        %                 acor(uu,ua) = abs(sys.mtH(uu,:)*sys.mtH(ua,:)'/(norm(sys.mtH(ua,:))*norm(sys.mtH(uu,:))))^2;
        %                 %acor(ua,uuu) = acor(uu,ua);
        %             end
        %         end
        %         acor
        
        
        
        if (~isempty(strmatch('MA',iRelayScheme)))
            [ vtS ] = MAsinr( pa, sys );
        elseif (~isempty(strmatch('MG',iRelayScheme)))
            [ vtS ] = MGsinr( pa, sys );
            %snr (worst case)
            % rotação vai constante dos filtros não altera o resultado MG
            % sem relay
            %             rxPowVtr = zeros(1,iUsers);
            %             intePowVtr = rxPowVtr;
            %             %rxPow2Vtr = rxPowVtr;
            %             %intePow2Vtr = rxPowVtr;
            %             for uu=1:iUsers
            %                 gu = sys.ueGrpVtr(uu);
            %                 for gg=1:pa.numOfGroups
            %                     if gg == gu
            %                         rxPowVtr(uu) = (abs(sys.mtH(uu,:)*sys.M{gg}).^2);
            %                         %rxPow2Vtr(uu) = (abs(sys.mtHr(uu,:)*sys.M{gg}*sys.dDr).^2);
            %                     else
            %                         intePowVtr(uu) = intePowVtr(uu) + (abs(sys.mtH(uu,:)*sys.M{gg}).^2);
            %                         %intePow2Vtr(uu) = intePow2Vtr(uu) + (abs(sys.mtHro(uu)*sys.dDr).^2 + 1);
            %
            %                         %(abs(sys.mtH(uu,:)*vtM{gg}).^2);
            %
            %                         %intePowVtr(uu)*pa.dBasePower
            %                     end
            %                 end
            %             end
            %             snrOne = (rxPowVtr*pa.dBasePower)./((intePowVtr*pa.dBasePower) + pa.dNoisePwr);
            %             %snrTwo = (rxPow2Vtr*pa.dBasePower)./((intePow2Vtr*pa.dBasePower) + pa.dNoisePwr);
            %             vtS = snrOne;%(snrOne + snrTwo);%*pa.dBasePower;
            %             %             dNpow = ((abs(sys.mtHro(uu)*dDr).^2 + 1)*pa.dNoisePwr);
            %             %             snrTwo(uu) = (abs(sys.mtHr(uu,:)*vtM*dDr).^2) / dNpow;
            
            %min(snrOne)
            %[abs(norm(sys.M{1}))^2  abs(norm(sys.M{2}))^2  abs(norm(sys.M{3}))^2]
        else
            %snr1;
            snrOne = snr1(sys,pa); %sys.mtHro,sys.M,sys.mtH
            %snr2;
            snrTwo = snr2(sys,pa);%sys.mtHro,sys.dDr,sys.M,sys.mtHr
            vtS = (snrOne + snrTwo)*pa.dBasePower;
        end
        
        vtSNRi = vtS;
        
        % store SNR for all iterations, sysmbols and users
        mtSNR(i,j,:) = vtSNRi;
        
    end
    
end

RandStream.setDefaultStream(mainStream);

stRes.mtSNR = mtSNR;

