function [AllSniffs, ColumnInfo, SingleUnits] = GetAllSniffs(MySession, varargin)

%% parse input arguments
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('quickmode', false, @(x) islogical(x) || x==0 || x==1);

% extract values from the inputParser
params.parse(varargin{:});
quickmode = params.Results.quickmode;

if ~quickmode
    load(MySession); % loads TracesOut PassiveOut SingleUnits
    max_sess = 2;
else
    TracesOut = MySession;
    max_sess = 1;
end

% get sniffs for both closed loop and passive phase
AllSniffs = [];
for sessionphase = 1:max_sess

    switch sessionphase
        case 1 % closed loop
            Traces = TracesOut;
        case 2
            Traces = PassiveOut;
    end

    %% Sniff Onsets and Offsets
    Sniffidx = find(diff(abs(Traces.SniffsDigitized{1})));
    Sniffidx = reshape(Sniffidx,2,[])';
    Sniffidx(:,1) = Sniffidx(:,1) +  1;
    Sniffidx(1:end-1,3) = Sniffidx(2:end,1); % next inhalation index
    Sniffidx(end,:) = [];

    %% separate various sniff types
    for n = 1:size(Sniffidx,1)

        % ITI or not
        ITIportion = ...
            numel(find( Traces.Manifold{1}(Sniffidx(n,1):Sniffidx(n,2)) )) / ...
            (Sniffidx(n,2) - Sniffidx(n,1) + 1 ); % fraction of time when manifold was on during inhalation
        % Odor identity
        StimType = mode( Traces.Odor{1}(Sniffidx(n,1):Sniffidx(n,2)) );
        Sniffidx(n,4) = ITIportion;
        Sniffidx(n,5) = StimType;
        Sniffidx(n,6) = mean(Traces.SniffsLocationed{1}(Sniffidx(n,1):Sniffidx(n,2))); % odor location
        Sniffidx(n,7) = mode( Traces.Trial{1}(Sniffidx(n,1):Sniffidx(n,2)) ); % perturbed or not

        Sniffidx(n,8) = sessionphase;

    end

    Sniffidx(:,11:13) = Sniffidx(:,1:3); % keep the indices
    Sniffidx(:,1:2) = Traces.Timestamps{1}(Sniffidx(:,1:2));
    Sniffidx(:,3) = Sniffidx(:,2) - Sniffidx(:,1);

    % give each sniff a event ID - which odor transition does it belong to,
    % and within that give a sniff rank
    OdorVector = Traces.Odor{1};
    OdorVector(OdorVector>0) = 1;
    OdorVector(OdorVector<1) = 0;
    OdorTransitions = [];
    OdorTransitions(:,1) = find(diff([0; OdorVector])>0); % first index when odor turns ON
    OdorTransitions(:,2) = find(diff([OdorVector; 0])<0); % last index when odor is still ON
    for Transition = 1:size(OdorTransitions,1)-1 % every trial
        % first sniff post odor start
        odorsniffs  = intersect(find(Sniffidx(:,11)>OdorTransitions(Transition,1)), ...
                               find(Sniffidx(:,11)<OdorTransitions(Transition,2)) );
        postsniffs  = intersect(find(Sniffidx(:,11)>OdorTransitions(Transition,2)), ...
                               find(Sniffidx(:,11)<OdorTransitions(Transition+1,1)) );
        % split Manifold on and Manifold off
        airsniffs = postsniffs(find( Sniffidx(postsniffs,4)),:);
        ITIsniffs = postsniffs(find(~Sniffidx(postsniffs,4)),:);

        Sniffidx(odorsniffs,9) = Transition;
        Sniffidx(odorsniffs,10) = 1;
        Sniffidx(airsniffs,9) = Transition;
        Sniffidx(airsniffs,10) = 0;
        Sniffidx(ITIsniffs,9) = Transition;
        Sniffidx(ITIsniffs,10) = -1;

    end


    AllSniffs = vertcat(AllSniffs, Sniffidx);

end

%SniffsOut.sniffs = AllSniffs;
ColumnInfo = {'InhStartTS', 'InhEndTS', 'InhDur', ...
                        'ManifoldState', 'OdorState', 'OdorLocation', ...
                        'PerturbationState', 'SessionPhase', ...
                        'event#', 'eventphase', ...
                        'InhStartIdx', 'InhEndIdx', 'NextInhIdx'};

end