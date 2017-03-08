function createOutputFiles(pa)

% para cada scheme
for s=1:length(pa.vtScheme),
    
    %strFileName = [pa.simulation '_' num2str(pa.sKstr) '_TX' num2str(pa.iTxAnt) '_UE' num2str(pa.iUsers) '_SC_' pa.vtScheme{s}];
    strFileName = par2filename(pa,s);
    
    fid = fopen([strFileName '.mat']);
    if fid == -1,
        %% se o arquivo não existe
        if ~isempty(strmatch('rotate',pa.vtScheme{s}))%strcmp(pa.vtScheme{s},'rotate') == 1,
            %% com 1 algoritmo
            if pa.getBER
                ou.mtBER = [];
            end
            ou.cellSNR = [];
            ou.cellAux = {};
            save(strFileName,'pa','ou');
        else
            %% com multiplos algoritmos
            if pa.getBER
                ou.mtBER = zeros(1,length(pa.vtAlgs));
            end
            ou.cellSNR = cell(1,length(pa.vtAlgs));
            ou.cellAux = {};           
            
            save(strFileName,'pa','ou');
        end
    else
        %% se o arquivo existe
        fclose(fid);
    end
end