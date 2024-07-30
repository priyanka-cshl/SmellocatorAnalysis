function [SniffTS] = ProcessThermistorData(RespirationData,varargin)

%% parse input arguments
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('samplerate', 500, @(x) isnumeric(x));
params.addParameter('plotting', 0, @(x) islogical(x));

% extract values from the inputParser
params.parse(varargin{:});
SampleRate = params.Results.samplerate;
plotting = params.Results.plotting;

% was odorlocation info provided?
if size(RespirationData,2) == 3
    OdorLocations = RespirationData(:,3);
    RespirationData(:,3) = [];
else
    OdorLocations = [];
end

% default settings for filtering
fband = [0.1 30];
Np    = 4; % filter order
[b,a] = butter(Np,fband/(SampleRate/2)); % band pass Butterworth filter coefficients

%%
% delaing with nans
invalid_indices = [];
if any(isinf(RespirationData(:,2)))
    x1 = find(isinf(RespirationData(:,2)),1,'first') - 1;
    x2 = find(isinf(RespirationData(:,2)),1,'last')  + 1;
    RespirationData(x1:x2,2) = linspace(RespirationData(x1,2),RespirationData(x2,2),(x2-x1+1));
    invalid_indices = x1:x2;
end

%% filtering - thermistor
TH_filt = filtfilt(b,a,RespirationData(:,2));
RespirationData(:,3) = TH_filt;

%% detect inflexion points
peakprom = std(TH_filt)/2.5;
% inhalation start?
[pks_ex,locs_ex] = findpeaks(TH_filt,'MinPeakProminence',peakprom);
% inhalation end
[pks_in,locs_in] = findpeaks(-TH_filt,'MinPeakProminence',peakprom);

% ignore any indices than were interpolated
if any(ismember(locs_ex,invalid_indices))
    whichones = find(ismember(locs_ex,invalid_indices));
    figure;
    temptrace = TH_filt;
    temptrace(invalid_indices) = NaN;
    plot(RespirationData(:,1),temptrace);
    hold on
    plot(RespirationData(locs_ex(whichones),1),pks_ex(whichones),'vk');
    todelete = userchoice(whichones);
    locs_ex(todelete,:) = [];
    pks_ex(todelete,:) = [];
    close(gcf);
end
if any(ismember(locs_in,invalid_indices))
    whichones = find(ismember(locs_in,invalid_indices));
    figure;
    temptrace = TH_filt;
    temptrace(invalid_indices) = NaN;
    plot(RespirationData(:,1),temptrace);
    hold on
    plot(RespirationData(locs_in(whichones),1),-pks_in(whichones),'.r');
    todelete = userchoice(whichones);
    locs_in(todelete,:) = [];
    pks_in(todelete,:) = [];
    close(gcf);
end

if plotting
    figure;
    plot(RespirationData(:,1),TH_filt);
    hold on
    plot(RespirationData(locs_ex,1),pks_ex,'vk');
    plot(RespirationData(locs_in,1),-pks_in,'.r');
end

SniffTS = [];
% sanity checks - missed peaks or double peaks
sortlocs = compute_sortlocs(locs_ex,locs_in);

while any(diff(sortlocs(:,3))==0)

    figure;
    plot(RespirationData(:,1),TH_filt);
    hold on
    plot(RespirationData(locs_ex,1),RespirationData(locs_ex,3),'.k');
    plot(RespirationData(locs_in,1),RespirationData(locs_in,3),'.r');
    
    f = find(diff(sortlocs(:,3))==0,1,'first');
    if sortlocs(f,2) > 0 % missing an inhalation detection
        % try local valley detection
        [~,missedval] = min(TH_filt(sortlocs(f,1):sortlocs(f+1,1)));
        missedidx = missedval - 1 + sortlocs(f,1);
        plot(RespirationData(missedidx,1),TH_filt(missedidx),'or');
        [locs_in,locs_ex] = editpoints(locs_in,missedidx,locs_ex);
    else % missing an exhalation detection
        % try local peak detection
        [~,missedpeak] = max(TH_filt(sortlocs(f,1):sortlocs(f+1,1)));
        missedidx = missedpeak - 1 + sortlocs(f,1);
        plot(RespirationData(missedidx,1),TH_filt(missedidx),'ok')
        [locs_ex,locs_in] = editpoints(locs_ex,missedidx,locs_in);
    end
    sortlocs = compute_sortlocs(locs_ex,locs_in);
    close(gcf);
end

firstinhalation = find(sortlocs(:,3)==0,1,"first");
for i = firstinhalation:2:size(sortlocs,1) % every inhalation
    if i > 1
        sniffstart = RespirationData(sortlocs(i-1,1),1);
    else
        sniffstart = NaN;
    end
    sniffend    = RespirationData(sortlocs(i,1),1);
    if i < size(sortlocs,1)
        nextsniff   = RespirationData(sortlocs(i+1,1),1);
    else
        nextsniff   = NaN;
    end
    if ~isempty(OdorLocations) && ~isnan(sniffstart)
        thisTrialOdorLocation = median(OdorLocations(sortlocs(i-1):sortlocs(i,1)));
        SniffTS     = vertcat(SniffTS, [sniffstart sniffend nextsniff thisTrialOdorLocation]);
    else
        SniffTS     = vertcat(SniffTS, [sniffstart sniffend nextsniff NaN]);
    end
end

%% function definitions
    function [todelete] = userchoice(whichones)
        answer = questdlg('Keep the identified point?', ...
            'Peak/Valley curation', ...
            'Yes','No','No');
        
        % Handle response
        switch answer
            case 'Yes'
                todelete = whichones;
            case 'No'
                todelete = [];
        end
    end

    function [locs_1,locs_2] = editpoints(locs_1,locs_to_edit,locs_2)
        
        list = { 'Keep new detection', 'Delete 1','Delete 2','Other'};
        [answer,tf] = listdlg('ListString',list);
        
        % Handle response
        if tf
            switch answer
                case 1
                    locs_1 = sort([locs_1; locs_to_edit]);
                case 2
                    locs_to_delete = find(locs_2<=locs_to_edit,1,'last');
                    locs_2(locs_to_delete,:) = [];
                case 3
                    locs_to_delete = find(locs_2>=locs_to_edit,1,'first');
                    locs_2(locs_to_delete,:) = [];
                case 4
                    keyboard;
            end
        end
    end

    function [sortlocs] = compute_sortlocs(locs_ex,locs_in)
        locs = vertcat([locs_ex(:,1) (1:size(locs_ex,1))'],[locs_in(:,1) -1*(1:size(locs_in,1))']);
        sortlocs = sortrows(locs,1);
        temp = sortlocs(:,2);
        temp(temp>0) = 1;
        temp(temp<0) = 0;
        sortlocs(:,3) = temp;
    end


end