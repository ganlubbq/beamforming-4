function simulador_sgb(pUE,pSimStr)
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%
% Correspondence concerning this software should be addressed to the e-mail
%brauner@virtual.ufc.br
%Ricardo Brauner dos Santos
%Professor
%Universidade Federal do Ceará
tic
addpath(genpath('.'));

pa = Parameters;  %cheat.pa();

%% read parametes
if nargin > 0
    pa.iUsers = pUE;
    pa.subsim = pSimStr;
end

%% parametros para variar
% parcela dos ues no hotspot
% potencia do relay

%% Definition of parameters to be varied

% NADA DA CODIFICAÇÃO POR LETRAS ABAIXO VALE AINDA!
%% the scheme name will be composed of prefixes that define:
% multicating scheme:
%% relay scheme
% SINGLE TRANSMISSION
% n - noRelay
% r - onlyRelay
% NO TRANSMIT PROCESSING
% e - eqCombine
% m - MRC
% INTELIGENT TRNASMISSION
% i - iterative
% s - scheme specific
%% number of antennas on relay
% b - basic, single antenna\
% a - relayed multi antenna multi group NOT IMPLEMENTED
%% number of groups
% s - single group
% m - multi group

%% micro algorithm
% 0, Matched filter
% 2, Zero forcing

% the fist letter must be S (caps) to describe an scheme that user a
% sub-algorithm

pa.vtScheme = {'test','MArand','SnoRelay'}; 
pa.vtAlgs = 0;%[8 10];%6 [0 6 8 10];
%[0 6 8 10 11];%[0 1 2 3 4 5 7 8];


disp(['XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX UEs: ' num2str(pa.iUsers) ' sub: ' pa.subsim]);

%% Sets the simulation state
vtSeed = getStates(pa);

%% Check whether saved files already exist, initializing them when necessary
createOutputFiles(pa)

toc
%% main loop
for s=1:length(pa.vtScheme),
    
    %stRes = cell(1,0);
    stRes = cell(1,length(pa.vtAlgs));
    stAux = cell(1,length(pa.vtAlgs));
    
    %% perform simulations
    %if (~isempty(strmatch('rotate',pa.vtScheme{s}))) || (~isempty(strmatch('MG',pa.vtScheme{s})))
    if (pa.vtScheme{s}(1)~='S')
        %% sets the coorect seed
        defaultStream = RandStream.getDefaultStream;
        defaultStream.State = vtSeed;
        
        %% displays simulations
        disp(['TxAnt ' num2str(pa.iTxAnt) ',scheme ' pa.vtScheme{s} ', No Alg.']);
        tic;
        %% runs the pontual simulation
        [stRes stAux] = linksim_sgb(pa, -1, pa.vtScheme{s});
        
        %% displays time
        toc;
        %stAux{1}
    else
        for i=1:length(pa.vtAlgs),
            
            %% sets the coorect seed
            defaultStream = RandStream.getDefaultStream;
            defaultStream.State = vtSeed;
            
            %% displays simulations
            disp(['TxAnt ' num2str(pa.iTxAnt) ',scheme ' pa.vtScheme{s} ', Alg ' num2str(pa.vtAlgs(i)) ' .']);
            tic;
            %% runs the pontual simulation
            [stRes{i} stAux{i}] = linksim_sgb(pa, pa.vtAlgs(i), pa.vtScheme{s});
            
            %% displays time
            toc;
            %stAux{1}
        end
    end
    
    %% loading data from previous runs and saving results
    %% TODO: compatibility now
    saveOutput(pa,stRes,stAux,s);
    
    disp('File saved ...');
end
toc
%figure;
%plot(1:5);
