function [SniffCoords,SniffProps] = TallyThermistorNPressureSniffs(WhereSession)

load(WhereSession,"SniffCoords","SniffProps");
if exist("SniffCoords")
    if ~isempty(SniffProps)
        return;
    end
else
    load(WhereSession, '-regexp','\w*SniffTimestamps$');
    SniffCoords = [];
    SniffProps = [];
end
if exist('CuratedMFSSniffTimestamps','var') & exist('CuratedSniffTimestamps','var')
    load(WhereSession,'RespirationData');

    mfs     = CuratedMFSSniffTimestamps;
    therm   = CuratedSniffTimestamps;
    % check for any missed thermistor peaks
    if any(mfs(:,7)<0)
        keyboard;
    end
    % are there any missed thermistor inhalations
    for i = 1:size(mfs,1)-1
        myTh = mfs(i,5);
        [deltaT,whichsniff] = min(abs(therm(:,1)-myTh));
        if deltaT>=0.003
            keyboard;
        end
        SniffCoords(i,1:5)  = therm(whichsniff,[1 2 3 8 9]); % [inhStart inhStop nextinhStart inhStartIdx inhStopIdx] - therm
        SniffProps(i,1:5)   = [i therm(whichsniff,[4 5 6 7])];
        SniffCoords(i,6:10) = mfs(i,[1 2 3 8 9]); % [inhStart inhStop nextinhStart inhStartIdx inhStopIdx] - predTherm
        
        % corresponding mfs zero crossing
        [deltaT,mfsidx] = min(abs(RespirationData(:,1)-mfs(i,4)));
        if deltaT>=0.003
            keyboard;
        end
        SniffCoords(i,11:15)= [RespirationData(mfsidx,1) mfs(i,2) mfs(i+1,4) mfsidx mfs(i,9)]; % [inhStart inhStop nextinhStart inhStartIdx inhStopIdx] - MFS zero crossings
    end

    % were any thermistor sniffs skipped
    if any(diff(SniffProps(:,1)~=-1))
        keyboard;
    end

    save(WhereSession,'SniffCoords','SniffProps','-append');
end
end