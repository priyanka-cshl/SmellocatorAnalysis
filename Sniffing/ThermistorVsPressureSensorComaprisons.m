% comparing latencies and sniff dynamics between 
% pressure sensor and thermistor signals

%% Load the data 
DataRoot = '/Users/Priyanka/Desktop/LABWORK_II/Data/Behavior/';
DataPath = 'Therm3/Therm3_20190927_r0.mat'; % also 28, 27
%DataPath = 'Therm2/Therm2_20190907_r0.mat'; % also 4,5,6,7,8,11 and 4 r1, r2

load(fullfile(DataRoot,DataPath));

PS = session_data.trace(:,5); % pressure sensor
TH = session_data.trace(:,6); % thermistor
TS = session_data.timestamps; % timestamps

%% Detect Inhalation exhalation timestamps   
% using the thermistor
[SniffTS_TH, inflexion_points_TH, THfilt] = ProcessThermistorData([TS TH]); %,'plotting',true);

% using the pressure sensor
[SniffTS_PS, inflexion_points_PS, PSfilt] = ProcessThermistorData([TS PS]); %,'plotting',true);

% getting pressure from thermistor
PS2 = diff(TH);
PS2(end+1) = PS2(end);
[SniffTS_PS2, inflexion_points_PS2, PSfilt2] = ProcessThermistorData([TS PS2]); %,'plotting',true);

%%
figure;
plot(TS,PSfilt,'k'); % pressure sensor
hold on;
plot(TS,-1.5+THfilt*10,'b'); % thermistor: offset and scaled
set(gca, 'XLim', [500 520]);

% peaks-valleys as detected from pressure sensor - plotted on PS trace
inflexion_points = inflexion_points_PS;
plot(TS(inflexion_points(find(~inflexion_points(:,3)),1)),...
    PSfilt(inflexion_points(find(~inflexion_points(:,3)),1)),'or');
plot(TS(inflexion_points(find(inflexion_points(:,3)),1)),...
    PSfilt(inflexion_points(find(inflexion_points(:,3)),1)),'ok');

% peaks-valleys as detected from pressure sensor - plotted on TH trace
plot(TS(inflexion_points(find(~inflexion_points(:,3)),1)),...
    -1.5 + 10*THfilt(inflexion_points(find(~inflexion_points(:,3)),1)),'.r');
plot(TS(inflexion_points(find(inflexion_points(:,3)),1)),...
    -1.5 + 10*THfilt(inflexion_points(find(inflexion_points(:,3)),1)),'.k');

% peaks-valleys as detected from thermistor - plotted on TH trace
inflexion_points = inflexion_points_TH;
plot(TS(inflexion_points(find(~inflexion_points(:,3)),1)),...
    -1.5 + 10*THfilt(inflexion_points(find(~inflexion_points(:,3)),1)),'or');
plot(TS(inflexion_points(find(inflexion_points(:,3)),1)),...
    -1.5 + 10*THfilt(inflexion_points(find(inflexion_points(:,3)),1)),'ok');

% peaks-valleys as detected from thermistor - plotted on PS trace
plot(TS(inflexion_points(find(~inflexion_points(:,3)),1)),...
    PSfilt(inflexion_points(find(~inflexion_points(:,3)),1)),'.r');
plot(TS(inflexion_points(find(inflexion_points(:,3)),1)),...
    PSfilt(inflexion_points(find(inflexion_points(:,3)),1)),'.k');

%% overlay pressure sensor from thermistor diff
figure;
plot(TS,-PSfilt2*150,'color',Plot_Colors('t'),'Linewidth',1);
hold on;
plot(TS,PSfilt,'k'); % pressure sensor
plot(TS,-1.5+THfilt*10,'b'); % thermistor: offset and scaled
set(gca, 'XLim', [500 520],'YLim',[-3 1.5]);

%% predicting PS from TH
figure;
scatter(PSfilt,-PSfilt2*150,'.k');
set(gca,'XLim',[-1.5 3.5],'YLim',[-1.5 3.5]);
axis square;

%% comparing the two sensors
Delta_starts = [];
for n = 1:size(SniffTS_TH,1)
    % find the valley just preceding each sniff start
    t_th = SniffTS_TH(n,1); % ts of sniff start in thermistor trace
    x1   = find(TS==t_th);
    t_ps = find(SniffTS_PS(:,2)<t_th,1,'last'); % valley closest to this peak in thermistor trace.
    x2   = find(TS==SniffTS_PS(t_ps,2));
    y1   = find(PSfilt(x2:end)>=0,1,'first');
    %line([TS(x1) TS(y1+x2-1)],[(-1.5 + 10*THfilt(x1)) PSfilt(y1+x2-1)],'LineStyle',':','color','r');
    Delta_starts(n,1) = TS(y1+x2-1) - TS(x1);
    Delta_starts(n,2) = SniffTS_TH(n,3) - SniffTS_TH(n,1); % sniff duration
    Delta_starts(n,3) = SniffTS_TH(n,2) - SniffTS_TH(n,1); % inhalation duration
end

figure;
subplot(2,2,1);
scatter(Delta_starts(:,2),Delta_starts(:,1),'.k');
subplot(2,2,2);
scatter(Delta_starts(:,3),Delta_starts(:,1),'.k');
subplot(2,2,3);
scatter(Delta_starts(:,2),Delta_starts(:,1),'.k');
set(gca,'YLim',[-0.05 0.15]);
subplot(2,2,4);
scatter(Delta_starts(:,3),Delta_starts(:,1),'.k');
set(gca,'YLim',[-0.05 0.15]);

%% for the thermistor: plotting sniff features vs. duration
SniffFeatures = [];
for n = 1:size(SniffTS_TH,1)
    x1 = find(TS==SniffTS_TH(n,1));
    x2 = find(TS==SniffTS_TH(n,2));
    delta_inhalation = SniffTS_TH(n,2) - SniffTS_TH(n,1);
    delta_sniff = SniffTS_TH(n,3) - SniffTS_TH(n,1);
    delta_amplitude = TH(x1) - TH(x2);
    inhalation_slope = delta_amplitude/delta_inhalation;
    SniffFeatures(n,:) = [delta_inhalation delta_sniff delta_amplitude inhalation_slope THfilt(x2)];
    %plot(delta_time,delta_depth,'ok');
end

figure; 
subplot(2,3,1); % plot inhalation duration vs inhalation peak
scatter(SniffFeatures(:,1),SniffFeatures(:,5),'.k');
subplot(2,3,2); % plot inhalation duration vs inhalation amplitude
scatter(SniffFeatures(:,1),SniffFeatures(:,3),'.k');
subplot(2,3,3); % plot inhalation duration vs inhalation slope
scatter(SniffFeatures(:,1),SniffFeatures(:,4),'.k');
subplot(2,3,4); % plot inhalation duration vs inhalation peak
scatter(SniffFeatures(:,2),SniffFeatures(:,5),'.k');
subplot(2,3,5); % plot inhalation duration vs inhalation amplitude
scatter(SniffFeatures(:,2),SniffFeatures(:,3),'.k');
subplot(2,3,6); % plot inhalation duration vs inhalation slope
scatter(SniffFeatures(:,2),SniffFeatures(:,4),'.k');

%% Plotting sniff waveforms 
% resort by inhalation duration
[~,sortorder] = sort(SniffTS_TH(:,3)-SniffTS_TH(:,1));
SniffsSorted = SniffTS_TH(sortorder,:);
nSniffs = size(SniffTS_TH,1);

SniffMat = ones(nSniffs,ceil(max(SniffsSorted(:,3)-SniffsSorted(:,1))/.002));
SniffOffs = zeros(nSniffs,1);

sniffcolors = brewermap(round(1.5*nSniffs),'*RdPu');

figure;
subplot(1,2,1);
hold on
for s = 1:nSniffs % every sniff
    x1 = find(TS==SniffsSorted(s,1));
    x2 = find(TS==SniffsSorted(s,2));
    x3 = find(TS==SniffsSorted(s,3));

    if ~isempty(x1) & ~isempty(x3)
        % plot the thermistor
        thissniff = THfilt(x1:x3);
        thissniff = thissniff - thissniff(1); % start all of them at zero amplitude
        %thissniff = thissniff/abs(min(thissniff)); % normalize to sniff amplitude
        plot3(TS(x1:x3)-TS(x1),0.1*s + 0*TS(x1:x3),thissniff,'color',sniffcolors(s+round(0.4*nSniffs),:));
        
        binindices = ceil((TS(x1:x3)-TS(x1))/0.002);
        binindices(1) = 1;
        SniffMat(s,binindices) = thissniff;
        SniffOffs(s,1) = numel(thissniff);
    end
end
view([-0.1 -0.4 0.2]);
set(gca,'XLim',[0 0.7],'YLim',[1 650],'ZLim',[-0.25 0.025]);

subplot(1,2,2);
imagesc(SniffMat(:,1:350),[-0.175 -0.0125]);
line(reshape([SniffOffs SniffOffs circshift(SniffOffs,-1)]',[],1),reshape([(1:nSniffs)'-1 (1:nSniffs)' nan*(1:nSniffs)']',[],1),'color','k','Linewidth',1);
set(gca,'Colormap',brewermap([],'*PuRd'));

%%
Blocks = round(linspace(1,nSniffs+1,4));
Blocks(2,1:end-1) = Blocks(1,2:end)-1;
Blocks(:,end) = [];

% for plotting
colors = brewermap(round(nSniffs/2),'*RdPu');

figure; 
for s = 1:numel(Blocks)
    subplot(2,size(Blocks,2),s);
    hold on;
end

for b = 1:size(Blocks,2)
    % every sniff in this duration block
    for n = Blocks(1,b):Blocks(2,b)
        x2 = find(TS==SniffsSorted(n,1));
        x2 = find(TS==SniffsSorted(n,2));
        x3 = find(TS==SniffsSorted(n,3));
        
        if ~isempty(x2) & ~isempty(x3)
            % plot the thermistor
            thissniff = TH(x2:x3);
            thissniff = thissniff - thissniff(1);
            %thissniff = thissniff/abs(min(thissniff));

            subplot(2,size(Blocks,2),b);
            plot(thissniff,'color',colors(n-Blocks(1,b)+1,:));

            % plot the pressure sensor
            thissniff = PS(x2:x3);
            thissniff = thissniff - thissniff(1);

            %thissniff = thissniff/abs(min(thissniff));
            subplot(2,size(Blocks,2),b+size(Blocks,2));
            plot(thissniff,'color',colors(n-Blocks(1,b)+1,:));
        end
    end
end

%%

% 
%     
% figure;
% plot(-PS/10);
% hold on
% plot(TH);
% 
% for n = 1:numel(locs_ex)
%     line([locs_ex(n) locs_ex(n)], [-0.2 0.2], 'Color', 'k', 'LineStyle', ':');
% end
% 
% for n = 1:numel(locs_in)
%     line([locs_in(n) locs_in(n)], [-0.2 0.2], 'Color', 'r', 'LineStyle', ':');
% end
% 
% set(gca,'XLim',10^5.*[2.95 2.97],'YLim',[-0.25 0.25]);
% set(gcf,'Renderer','painters');
% print('/Users/Priyanka/Desktop/LABWORK_II/Presentations/InHouse2024/ThermistorPressureSensorExample.eps','-depsc','-tiff','-r300','-painters');
% 
