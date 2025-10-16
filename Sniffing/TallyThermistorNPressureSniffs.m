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
%     % check for any missed thermistor peaks that were hand picked during
%     % curation
%     if any(mfs(:,7)<0)
%         keyboard;
%     end
    % are there any missed thermistor inhalations
    for i = 1:size(mfs,1)-1
        if mfs(i,7) ~= -1 && ~isnan(mfs(i,5))
           myTh = mfs(i,5); % time-point when a peak occured in the thermistor side
           [deltaT,whichsniff] = min(abs(therm(:,1)-myTh)); %  find the row entry in thermistor sniffs that correspond to the matched sniff
            if deltaT>=0.003
                keyboard;
            end
            SniffCoords(i,1:5)  = therm(whichsniff,[1 2 3 8 9]); % [inhStart inhStop nextinhStart inhStartIdx inhStopIdx] - therm
            SniffProps(i,1:5)   = [i therm(whichsniff,[4 5 6 7])];
        else
            rescued = 0;
            if ~isnan(mfs(i+1,5)) && ~isnan(mfs(i-1,5)) % ignored but can be assigned
                whichsniff = find((therm(:,1)>mfs(i-1,5))&(therm(:,1)<mfs(i+1,5)),1,"first");
                if ~isempty(whichsniff)
                    SniffCoords(i,1:5)  = therm(whichsniff,[1 2 3 8 9]);
                    SniffProps(i,1:5)   = [-i therm(whichsniff,[4 5 6 7])];
                    rescued = 1;
                end
            end
            if ~rescued 
                if mfs(i,7) == -1
                    inhduration = mfs(i,2) - mfs(i,1);
                    nextinh = find(mfs(:,5)>(mfs(i,5)+inhduration),1,"first");
                    SniffCoords(i,1:3)  = [mfs(i,5) mfs(i,5)+inhduration mfs(nextinh,5)];
                    [~,SniffCoords(i,4)] = min(abs(RespirationData(:,1)-SniffCoords(i,1)));
                    [~,SniffCoords(i,5)] = min(abs(RespirationData(:,1)-SniffCoords(i,1)));
                    SniffProps(i,1:5) = [-i nan(1,4)];
                else
                    %keyboard;
                    SniffCoords(i,1:5)  = nan(1,5); % [inhStart inhStop nextinhStart inhStartIdx inhStopIdx] - therm
                    SniffProps(i,1:5)   = nan(1,5);
                end
            end
        end
        SniffCoords(i,6:10) = mfs(i,[1 2 3 8 9]); % [inhStart inhStop nextinhStart inhStartIdx inhStopIdx] - predTherm
        % corresponding mfs zero crossing
        [deltaT,mfsidx] = min(abs(RespirationData(:,1)-mfs(i,4)));
        if deltaT>=0.003
            keyboard;
        end
        SniffCoords(i,11:15)= [RespirationData(mfsidx,1) mfs(i,2) mfs(i+1,4) mfsidx mfs(i,9)]; % [inhStart inhStop nextinhStart inhStartIdx inhStopIdx] - MFS zero crossings
    end

    badstretches = find(isnan(SniffCoords(:,4))|isnan(SniffCoords(:,9))|isnan(SniffCoords(:,14)));
    SniffCoords(badstretches,[1:3 6:8 11:13]) = -abs(SniffCoords(badstretches,[1:3 6:8 11:13]));

    % were any thermistor sniffs skipped
%     if any(diff(SniffProps(:,1)~=-1))
%         keyboard;
%     end

    save(WhereSession,'SniffCoords','SniffProps','-append');
end
end