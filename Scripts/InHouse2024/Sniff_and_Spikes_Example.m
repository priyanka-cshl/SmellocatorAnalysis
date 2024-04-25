%% script to plot Spikes and Respiration and trial times
global SampleRate

SetupSmellocatorGlobals;

% load the data
DataRoot = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/Processed/Behavior/';

DataPath = 'S12/S12_20230731_r0_processed.mat'; 
MyUnits = [58];
sniffoffset = 0;
BehaviorSession = '/Users/Priyanka/Desktop/LABWORK_II/Data/Behavior/S12/S12_20230731_r0.mat';

% DataPath = 'Q3/Q3_20221019_r0_processed.mat'; 
% MyUnits = [16 18];

% DataPath = 'O3/O3_20210927_r0_processed.mat'; 
% MyUnits = [2 3 7 9 11 13 14 15 19 25 27 28 33 35 42 43 45 48 54 58 59 62 66 72];

% DataPath = 'Q9/Q9_20221116_r0_processed.mat'; 
% MyUnits = [1 11 15 18 19 23 28 29 36 39 43 49 94]; % Q9

% DataPath = 'Q8/Q8_20221204_r0_processed.mat'; 
% MyUnits = [49]; % Q9
% sniffoffset = 0;

%% plot sniffing 
Temp = load(BehaviorSession,'session_data');
whichcol = find(ismember(Temp.session_data.trace_legend,'thermistor'));
TS = Temp.session_data.timestamps;    
% filtering - thermistor
TH = Temp.session_data.trace(:,whichcol);
fband = [0.1 30];
Np    = 4; % filter order
[b,a] = butter(Np,fband/(SampleRate/2)); % band pass Butterworth filter coefficients
TH_filt = filtfilt(b,a,TH);
    
% detect inflexion points
peakprom = std(TH_filt)/2;
% inhalation start?
[pks_ex,locs_ex] = findpeaks(TH_filt,'MinPeakProminence',peakprom);
% inhalation end
[pks_in,locs_in] = findpeaks(-TH_filt,'MinPeakProminence',peakprom);

figure;
plot(TS,TH_filt); % filtered thermistor data
hold on
% plot(RespirationData(locs_ex,1),pks_ex,'vk');
% plot(RespirationData(locs_in,1),-pks_in,'.r');

for n = 1:numel(locs_ex)
    line([TS(locs_ex(n)) TS(locs_ex(n))], [-0.4 0.4], 'Color', 'b', 'LineStyle', ':');
end

for n = 1:numel(locs_in)
    line([TS(locs_in(n)) TS(locs_in(n))], [-0.4 0.4], 'Color', 'r', 'LineStyle', ':');
end

% plot trial times
whichcol = find(ismember(Temp.session_data.trace_legend,'trial_on'));
Trial = Temp.session_data.trace(:,whichcol)/10;
plot(TS,Trial);

% load units
load(fullfile(DataRoot,DataPath), 'TimestampAdjust', 'SingleUnits', 'SniffTS');

%% plot spiketimes of selected units
for n = 1:numel(MyUnits)
    spikes = SingleUnits(MyUnits(n)).spikes - TimestampAdjust.ClosedLoop; % in behavior time base
    PlotRaster(spikes,n-1,'k',0.3);
end

set(gca,'XLim',[1315 1322]);
set(gcf,'Renderer','painters');
print('/Users/Priyanka/Desktop/LABWORK_II/Presentations/InHouse2024/FrequencyOrOdorExample_1_lines.eps','-depsc','-tiff','-r300','-painters');

set(gca,'XLim',[1049 1056]);
set(gcf,'Renderer','painters');
print('/Users/Priyanka/Desktop/LABWORK_II/Presentations/InHouse2024/FrequencyOrOdorExample_2_lines.eps','-depsc','-tiff','-r300','-painters');
