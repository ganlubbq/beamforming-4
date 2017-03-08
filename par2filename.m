function name = par2filename(pa,iScheme)

name = [pa.simulation '_TX' num2str(pa.iTxAnt) '_UE' num2str(pa.iUsers) '_SC_' pa.vtScheme{iScheme}];