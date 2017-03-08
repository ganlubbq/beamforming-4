function [vtSeed] = getStates(pa)
%function [statesCell] = getStates(pa)

strFileName = ['states_' pa.simulation '_' pa.subsim '.mat'];
fid = fopen(strFileName);
if fid == -1,
%     % initializes to  first random seed
%     RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
%     % generates one random seed to each iteration
%     vtSeed = rand(1,pa.iIter)*1e5;
    
    % initializes to  first random seed
    defaultStream = RandStream('mt19937ar','seed',sum(100*clock));
    RandStream.setDefaultStream(defaultStream);
%     for ii=1:pa.iIter
%         defaultStream = RandStream('mt19937ar','seed',vtSeed(ii));
%         
%         statesCell{ii} = defaultStream.State;
%     end
    vtSeed = defaultStream.State;
    %save(strFileName,'statesCell');
    save(strFileName,'vtSeed');
else
    fclose(fid);
    load(strFileName,'vtSeed')
    %defaultStream = RandStream('mt19937ar','seed',sum(100*clock));
    %defaultStream.State = vtSeed;
    %RandStream.setDefaultStream(defaultStream);
end