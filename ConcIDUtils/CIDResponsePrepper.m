function [StimSettings, TTLs, SingleUnits, AllSpikes, TrialWiseSniffs, SniffsPlot] = CIDResponsePrepper(myKsDir,varargin)
%% Input
%myKsDir = '/mnt/storage/Sorted/Q8/2022-12-19_15-41-55';
%myKsDir = '/media/priyanka/ABC-ntfs/EphysSorted/Q9/2022-12-15_16-28-22';

%% parse input arguments
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('allUnits', 0, @(x) isnumeric(x)); % -1: non-good, 0: good (default), 1: all
params.addParameter('minRate', 0, @(x) isnumeric(x)); % in Hz
params.addParameter('align2sniffs', 1, @(x) isnumeric(x)); % 1 = first sniff, 2 = second sniff
params.addParameter('whichSniffSensor', 'mfs', @(x) ischar(x)); % 'mfs', 'therm'

% extract values from the inputParser
params.parse(varargin{:});
UnitFilter = params.Results.allUnits; % -1: non-good, 0: good (default), 1: all
minRate = params.Results.minRate; % in Hz
% align to sniffs or not, and which sensor to use
align2sniffs = params.Results.align2sniffs; % 1 = first sniff, 2 = second sniff
whichSniffSensor = 1+ strcmp(params.Results.whichSniffSensor,'mfs'); % 1 = thermistor, 2 = MFS


%% Input
%myKsDir = '/mnt/storage/Sorted/Q8/2022-12-19_15-41-55';
%myKsDir = '/media/priyanka/ABC-ntfs/EphysSorted/Q9/2022-12-15_16-28-22';

%% some default Settings
% for new conc. series expts, consider both pulses
bothPulses = 1;
% for conc. series, group together concentrations of a given odor
stackConcs = 1;

%% load the data
clear KS4Units;
load(fullfile(myKsDir,'quickprocesssniffs.mat')); % sniff times, KS4Units
load(fullfile(myKsDir,'quickprocessOdorTTLs.mat'),'TTLs','StimSettings'); % odorTTLs
SingleUnits = LoadKS4Units(myKsDir,'minSpikes',minRate,'allUnits',UnitFilter);
%disp(['found ',num2str(size(SingleUnits,2)),' good, >0.25Hz units']);
nUnits = size(SingleUnits,2);

% for new CID experiments
if size(TTLs.Odor,1) == 2*size(TTLs.Trial,1)
    StimSettings.SessionType = 'newCID';
    if bothPulses
        TTLs.Trial(:,2) = TTLs.Trial(:,2) + 2;
    end
end

%% stimulus settings for plotting
if stackConcs == 2
    if strcmp(StimSettings.SessionType,'ConcentrationSeries')
        newStimCol = (TTLs.Trial(:,5)*10^4) + TTLs.Trial(:,4);
        newTypesCol = TTLs.Trial(:,5)*0;
        TTLs.Trial(:,4:5) = [newStimCol newTypesCol];
        StimSettings.SessionType = '16_Concs*';
    end
end
if stackConcs == 1
    if strcmp(StimSettings.SessionType,'ConcentrationSeries')
        newStimCol = TTLs.Trial(:,5) + TTLs.Trial(:,4);
        newTypesCol = TTLs.Trial(:,5)*0;
        TTLs.Trial(:,4:5) = [newStimCol newTypesCol];
        StimSettings.SessionType = '16_Concs';
    end
end

nStim = unique(TTLs.Trial(:,4)); % no. of odors delivered
nTypes = unique(TTLs.Trial(:,5)); % concentrations used
StimSettings.nStim = nStim;
StimSettings.nTypes = nTypes;
StimSettings.align2sniffs = align2sniffs;

%% Hacks for some older experiments
% make air trials as last stimulus
TTLs.Trial((TTLs.Trial(:,4)==0),4) =  max(TTLs.Trial(:,4)) + 1;

if any(nStim == -1)
    f = find(TTLs.Trial(:,4)==-1);
    TTLs.Trial(f,5) = nTypes(2);
    nTypes(1,:) = [];
    nStim(nStim<0,:) = [];
    TTLs.Trial(f,7) = TTLs.Trial(f,1) + StimSettings.timing(2)/1000;
    TTLs.Trial(f,8) = TTLs.Trial(f,7) + StimSettings.timing(3)/1000;
end

if ~isfield(StimSettings,'SessionType')
    keyboard;
end

%% process sniffs
TrialWiseSniffs = [];
if exist('CuratedSniffTimestamps','var') || exist('CuratedMFSSniffTimestamps','var')
    addSniffPlot = 1;
end
if align2sniffs || addSniffPlot
    % find first inhalation start after odor start
    if whichSniffSensor==1 && exist('CuratedSniffTimestamps','var')
        MySniffTimeStamps = CuratedSniffTimestamps(:,1:3);
    elseif whichSniffSensor==2 && exist('CuratedMFSSniffTimestamps','var')
        MySniffTimeStamps = CuratedMFSSniffTimestamps(:,1:3);
    else
        keyboard;
    end
    for t = 1:size(TTLs.Trial,1)
        if TTLs.Trial(t,4)>0 % every valid trial
            ts = TTLs.Trial(t,[1 2 7 8]); % trial start, stop, odor start, stop
            if ~strcmp(StimSettings.SessionType,'newCID')
                if t < size(TTLs.Trial,1)
                    ts(2) = TTLs.Trial(t+1,1);
                else
                    ts(2) = ts(2) + (StimSettings.timing(5)/1000);
                end
            else
                disp('decide what to do here');
                keyboard;
            end
            % find all sniffs
            thisTrialSniffs = intersect(find(MySniffTimeStamps(:,1)>=ts(1)),find(MySniffTimeStamps(:,1)<ts(2)));
            firstsniff = find(MySniffTimeStamps(thisTrialSniffs,1)>=ts(3),1,'first');
            firstpostsniff = find(MySniffTimeStamps(thisTrialSniffs,1)>=ts(4),1,'first');
            TTLs.Trial(t,11) = MySniffTimeStamps(thisTrialSniffs(firstsniff+(align2sniffs-1)),1); % time of first sniff (true time)
            TTLs.Trial(t,12) = TTLs.Trial(t,7) - TTLs.Trial(t,11); % % time of first sniff from odor Onset
            thisTrialSniffs = MySniffTimeStamps(thisTrialSniffs,1:3) - ts(3); % odor start
            % now lets index every sniff
            thisTrialSniffs(:,4) = t;
            aftersniffs = find(thisTrialSniffs(:,1)>=0);
            beforesniffs = 1:(aftersniffs(1)-1);
            thisTrialSniffs(beforesniffs,5) = -flipud(beforesniffs(:));
            thisTrialSniffs(aftersniffs,5)  = (1:numel(aftersniffs))-1;
            thisTrialSniffs(:,6) = ts(3); % useful to get back actual value
            
            % add sniff phase
            thisTrialSniffs(beforesniffs,7) = 0;
            thisTrialSniffs(aftersniffs,7) = 1;
            postodorsniffs = find(thisTrialSniffs(:,1)>=(ts(4)-ts(3)));
            thisTrialSniffs(postodorsniffs,7) = 1.5;
            if strcmp(StimSettings.SessionType,'newCID')
                secondOdor = TTLs.Trial(t,[9 10]) - ts(3);
                postOdorTwo = find(thisTrialSniffs(:,1)>=(secondOdor(1)));
                thisTrialSniffs(postOdorTwo,7) = 2;
                postOdorTwo = find(thisTrialSniffs(:,1)>=(secondOdor(2)));
                thisTrialSniffs(postOdorTwo,7) = 2.5;
            end
            TrialWiseSniffs = vertcat(TrialWiseSniffs, thisTrialSniffs);
        end
    end
end

% TTLs.Trial has now sniff info in additional columns
% 11    : true time of first (or nth) sniff after odor ON
% 12    : relative time of first (or nth) sniff after odor ON

% TrialWiseSniffs: Columns
% 1 to 3: inh start, end, next w.r.t. this Trial's odor Onset
% 4     : trial index
% 5     : sniff index within a trial, 0 = first sniff after odor onset
% 6     : actual odor ON time
% 7     : trial phase : 0 - pre-odor, 1 - odor, 2 - second odor pulse,
%                       1.5 - post-odor, 2.5 - post second odor pulse
%                       3 - ITI

%% Make  trial-aligned Spike Plot
for n = 1:nUnits
    SpikesPlot = [];
    thisUnitSpikes = SingleUnits(n).spikes;
    for t = 1:size(TTLs.Trial,1)
        if TTLs.Trial(t,4)>0
            % every trial
            ts = TTLs.Trial(t,[1 2 7 8]); % trial start, stop, odor start, stop
            if size(TTLs.Trial,2) > 10 && align2sniffs
                ts(5) = TTLs.Trial(t, 11);
            else
                ts(5) = ts(3);
            end
            thistrialspikes = intersect(find(thisUnitSpikes>=ts(1)),find(thisUnitSpikes<ts(2)));
            thistrialspikes = thisUnitSpikes(thistrialspikes) - ts(5);
            multiplier = t;
            SpikesPlot = vertcat(SpikesPlot, ...
                [thistrialspikes multiplier*ones(numel(thistrialspikes),1)]);
        end
    end
    AllSpikes(n) = {SpikesPlot};
end

if addSniffPlot
    % add a fake unit with sniff times instead of spikes
    SpikesPlot = [];
    thisUnitSpikes = TrialWiseSniffs(:,1) + TrialWiseSniffs(:,6);
    for t = 1:size(TTLs.Trial,1)
        if TTLs.Trial(t,4)>0
            % every trial
            ts = TTLs.Trial(t,[1 2 7 8]); % trial start, stop, odor start, stop
            if size(TTLs.Trial,2) > 10 && align2sniffs
                ts(5) = TTLs.Trial(t, 11);
            else
                ts(5) = ts(3);
            end
            thistrialspikes = intersect(find(thisUnitSpikes>=ts(1)),find(thisUnitSpikes<ts(2)));
            thistrialspikes = thisUnitSpikes(thistrialspikes) - ts(5);
            multiplier = t;
            SpikesPlot = vertcat(SpikesPlot, ...
                [thistrialspikes multiplier*ones(numel(thistrialspikes),1)]);
        end
    end
    SniffsPlot = SpikesPlot;
else
    SniffsPlot = [];
end


