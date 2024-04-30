% comparing pressure sensor and thermistor signals

% filters
SampleRate = 500;
fband = [0.1 30];
Np    = 4; % filter order
[b,a] = butter(Np,fband/(SampleRate/2)); % band pass Butterworth filter coefficients

% load the data
DataRoot = '/Users/Priyanka/Desktop/LABWORK_II/Data/Behavior/';
DataPath = 'Therm3/Therm3_20190927_r0.mat'; % also 28, 27
%DataPath = 'Therm2/Therm2_20190907_r0.mat'; % also 4,5,6,7,8,11 and 4 r1, r2

load(fullfile(DataRoot,DataPath));

PS = session_data.trace(:,5); % pressure sensor
TH = session_data.trace(:,6); % thermistor

PS = filtfilt(b,a,PS);
TH = filtfilt(b,a,TH);

%% detect inflexion points
peakprom = std(TH)/2;
% inhalation start?
[pks_ex,locs_ex] = findpeaks(TH,'MinPeakProminence',peakprom);
% inhalation end
[pks_in,locs_in] = findpeaks(-TH,'MinPeakProminence',peakprom);
    
figure;
plot(-PS/10);
hold on
plot(TH);

for n = 1:numel(locs_ex)
    line([locs_ex(n) locs_ex(n)], [-0.2 0.2], 'Color', 'k', 'LineStyle', ':');
end

for n = 1:numel(locs_in)
    line([locs_in(n) locs_in(n)], [-0.2 0.2], 'Color', 'r', 'LineStyle', ':');
end

set(gca,'XLim',10^5.*[2.95 2.97],'YLim',[-0.25 0.25]);
set(gcf,'Renderer','painters');
print('/Users/Priyanka/Desktop/LABWORK_II/Presentations/InHouse2024/ThermistorPressureSensorExample.eps','-depsc','-tiff','-r300','-painters');

