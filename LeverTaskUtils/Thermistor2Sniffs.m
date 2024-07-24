function [SniffTS] = Thermistor2Sniffs(InputData,plotting)

if nargin<2
    plotting = 0;
end

RespirationData(:,1:2)  = InputData(:,1:2);
LocationData(:,1)       = InputData(:,3);

%% filtering - thermistor
TH = RespirationData(:,2);
SampleRate = 500;
fband = [0.1 30];
Np    = 4; % filter order
[b,a] = butter(Np,fband/(SampleRate/2)); % band pass Butterworth filter coefficients
TH_filt = filtfilt(b,a,TH);
RespirationData(:,3) = TH_filt;

%% detect inflexion points
peakprom = std(TH_filt)/2;
% inhalation start?
[pks_ex,locs_ex] = findpeaks(TH_filt,'MinPeakProminence',peakprom);
% inhalation end
[pks_in,locs_in] = findpeaks(-TH_filt,'MinPeakProminence',peakprom);

%     % detect exhalation start?
%     fband = [0.1 20];
%     [b,a] = butter(Np,fband/(SampleRate/2)); % band pass Butterworth filter coefficients
%     d_TH_filt = filtfilt(b,a,diff(TH_filt));
%     dd_TH_filt = diff(d_TH_filt);
%     dd_TH_filt = [dd_TH_filt; dd_TH_filt(end-1:end)];
%     % peaks in the first - derivative - possibly exhalation start and
%     % inhalation ends
%     [~ , locs_d] = findpeaks(d_TH_filt,'MinPeakProminence',std(d_TH_filt(d_TH_filt>0))/2);
%     [~ , locs_dd] = findpeaks(dd_TH_filt,'MinPeakProminence',std(dd_TH_filt)/4);

if plotting
    figure;
    plot(RespirationData(:,1),TH_filt);
    hold on
    plot(RespirationData(locs_ex,1),pks_ex,'vk');
    plot(RespirationData(locs_in,1),-pks_in,'.r');
end

SniffTS = [];
%%
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

% make the list
FirstExhalation = find(sortlocs(:,3)==1,1,'first');
LastExhalation = find(sortlocs(:,3)==1,1,'last') - 2;
for i = FirstExhalation:2:LastExhalation
    odorlocation = median(LocationData(sortlocs(i,1):sortlocs(i+1,1),1));
    SniffTS      = vertcat(SniffTS, ...
                   [RespirationData(sortlocs(i:i+2),1)' odorlocation] );
end

%%
% sanity checks - missed peaks or double peaks
% for i = 1:numel(locs_in)-2
%     putatives = intersect(find(locs_ex>=locs_in(i)), find(locs_ex<=locs_in(i+1)));
%     if numel(putatives) > 1
%         % multiple detections
%         figure;
%         plot(RespirationData((locs_in(i)-SampleRate):(locs_in(i+1)+SampleRate),2));
%         hold on
%         plot(RespirationData(locs_in(i),1),-pks_in(i),'.r');
%         plot(RespirationData(locs_in(i+1),1),-pks_in(i+1),'.r');
%         plot(RespirationData(locs_ex(putatives(1)),1),TH_filt(locs_ex(putatives(1)),1),'vk');
%         plot(RespirationData(locs_ex(putatives(2)),1),TH_filt(locs_ex(putatives(2)),1),'vk');
%         keyboard;
%     elseif numel(putatives) == 0
%         % missed points
%         figure;
%         plot(RespirationData((locs_in(i)-SampleRate):(locs_in(i+1)+SampleRate),2));
%         hold on
%         plot(RespirationData(locs_in(i),1),-pks_in(i),'.r');
%         plot(RespirationData(locs_in(i+1),1),-pks_in(i+1),'.r');
%         keyboard;
%     else
%         odorlocation = median(LocationData(locs_ex(putatives(1)):locs_in(i+1),1));
%         SniffTS = vertcat(SniffTS, [RespirationData(locs_ex(putatives(1)),1) RespirationData(locs_in(i+1),1) RespirationData(locs_ex(putatives(1)+1),1) odorlocation]);
%         % [inhalation-start inhalation-end next-inhalation]
%     end
% end

%%
if plotting
    % some pretty plots
    figure;
    [~,sortorder] = sort(SniffTS(:,3)-SniffTS(:,1));
    SniffTS_ = SniffTS(sortorder,:);
    line(repmat(1:size(SniffTS_,1),2,1), [(SniffTS_(:,1) - SniffTS_(:,1))'; (SniffTS_(:,3) - SniffTS_(:,1))'], 'color' ,'k');
    hold on
    line(repmat(1:size(SniffTS_,1),2,1), [(SniffTS_(:,1) - SniffTS_(:,1))'; (SniffTS_(:,2) - SniffTS_(:,1))'], 'color' ,'r');
    
    figure;
    SniffTS_ = SniffTS;
    line(repmat(1:size(SniffTS_,1),2,1), [(SniffTS_(:,1) - SniffTS_(:,1))'; (SniffTS_(:,3) - SniffTS_(:,1))'], 'color' ,'k');
    hold on
    line(repmat(1:size(SniffTS_,1),2,1), [(SniffTS_(:,1) - SniffTS_(:,1))'; (SniffTS_(:,2) - SniffTS_(:,1))'], 'color' ,'r');
    
    figure;
    figure, scatter(SniffTS(:,3) - SniffTS(:,1),SniffTS(:,2) - SniffTS(:,1), 'ok');
    set(gca,'YLim', [0 max(get(gca,'XLim'))]);
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