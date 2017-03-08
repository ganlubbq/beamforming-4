function saveOutput(pa,stRes,stAux,iScheme)

%% load previous output data
%strFileName = [pa.simulation '_' num2str(pa.sKstr) '_TX' num2str(pa.iTxAnt) '_UE' num2str(pa.iUsers) '_SC_' pa.vtScheme{iScheme}];
strFileName = par2filename(pa,iScheme);
loa = load(strFileName);

ou.mtBER = zeros(size(stRes));
ou.cellSNR = {};
ou.cellAux = {};
%ou.mtIter = ou.mtBER;

%% process output data
if ~iscell(stRes)%(~isempty(strmatch('rotate',pa.vtScheme{iScheme}))) || (~isempty(strmatch('MG',pa.vtScheme{iScheme})))
    %~isempty(strmatch('rotate',pa.vtScheme{iScheme})) %strcmp(,'rotate') == 1,
    ou.cellSNR = mean(min(stRes.mtSNR,[],3),2);
    
    ou.cellSNR = [ou.cellSNR;loa.ou.cellSNR];   
    
    if ~isempty(stAux);
        theAux = stAux{1};
    else
        theAux = stAux;
    end
    
    if iscell(theAux)
        if isempty(loa.ou.cellAux)
            
            for au=1:length(theAux)
                ou.cellAux{1}{au} = theAux{au};
            end
            
        else
            
            for au=1:length(theAux)
                ou.cellAux{1}{au} = [theAux{au};loa.ou.cellAux{1}{au}];
            end
        end
    else
        disp('theAux is not cell!');
        keyboard;
    end
else
    for alg=1:length(stRes)
        if pa.getBER
            ou.mtBER(alg) = mean(stRes{alg}.vtBER);
        end
        %% AKI
        ou.cellSNR{alg} = mean(min(stRes{alg}.mtSNR,[],3),2);
        %ou.mtIter(alg) = pa.iIter;
    end
    
    % store SNR results
    for alg=1:length(stRes)
        ou.cellSNR{alg} = [ou.cellSNR{alg};loa.ou.cellSNR{alg}];
        
        theAux = stAux{alg};
        
        if iscell(theAux)
            if isempty(loa.ou.cellAux)
                
                for au=1:length(theAux)
                    ou.cellAux{alg}{au} = theAux{au};
                end
                
            else
                
                for au=1:length(theAux)
                    ou.cellAux{alg}{au} = [theAux{au};loa.ou.cellAux{alg}{au}];
                end
            end
        else
            disp('theAux is not cell!');
            keyboard;
        end
    end
end

%% save output data
%mtTotalIter = loa.mtTotalIter + mtIter;

if pa.getBER
    % store BER results
    %re.mtBER = (loa.mtBER .* loa.mtTotalIter  + re.mtBER .* mtIter) ./ mtTotalIter;
end



if isfield(loa, 'sets')
    sets = loa.sets;
    sets{end+1} = pa.subsim;
else
    sets = {pa.subsim};
end

save(strFileName,'pa','ou','sets');