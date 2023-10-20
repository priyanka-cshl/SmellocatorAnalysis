function [PhasesOut] = WhichSniffPhase(EventsIn,SniffsIn,varargin)

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('warpmethod', 3, @(x) isnumeric(x)); % 0 - don't align, 1 - simple warp, 2 - complex warp, 3 - fixed latency

params.parse(varargin{:});
warpmethod = params.Results.warpmethod;

PhasesOut = EventsIn*0;
for i = 1:numel(EventsIn)
    whichsniff = find(SniffsIn(:,1)<=EventsIn(i),1,'last');
    if ~isempty(whichsniff)
        switch warpmethod % assymetric
            case 1
                thiseventphase = (EventsIn(i) - SniffsIn(whichsniff,1)) / ...
                    (SniffsIn(whichsniff,3) - SniffsIn(whichsniff,1));
                % from 0 to 1
            case 2
                % during inhalation or exhalation
                if EventsIn(i)<=SniffsIn(whichsniff,2)
                    % inhalation
                    thiseventphase = (EventsIn(i) - SniffsIn(whichsniff,1)) / ...
                        (SniffsIn(whichsniff,2) - SniffsIn(whichsniff,1));
                    thiseventphase = thiseventphase/2; % from 0 to 0.5
                else
                    % exhalation
                    thiseventphase = (EventsIn(i) - SniffsIn(whichsniff,2)) / ...
                        (SniffsIn(whichsniff,3) - SniffsIn(whichsniff,2));
                    thiseventphase = 0.5 + thiseventphase/2; % from 0.5 to 1
                end
            case 3
                thiseventphase = (EventsIn(i) - SniffsIn(whichsniff,1)) ; 
                % from 0 to 1 (in seconds)
        end
        PhasesOut(i) = SniffsIn(whichsniff,4) + thiseventphase;
    else
        PhasesOut(i) = NaN;
    end
end

end