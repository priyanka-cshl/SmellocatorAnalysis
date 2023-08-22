function [SniffTS] = ReadThermistorData(FileName)

%% load the file
Temp = load(FileName,'session_data');

if find(ismember(Temp.session_data.trace_legend,'thermistor'))
    whichcol = find(ismember(Temp.session_data.trace_legend,'thermistor'));
    RespirationData(:,2) = Temp.session_data.trace(:,whichcol);
    RespirationData(:,1) = Temp.session_data.timestamps;
    
    
    %% filtering - thermistor
    TH = RespirationData(:,2);
    SampleRate = 500;
    fband = [0.1 30];
    Np    = 4; % filter order
    [b,a] = butter(Np,fband/(SampleRate/2)); % band pass Butterworth filter coefficients
    TH_filt = filtfilt(b,a,TH);
    
    %% detect inflexion points
    peakprom = std(TH_filt)/2;
    % inhalation start?
    [pks_ex,locs_ex] = findpeaks(TH_filt,'MinPeakProminence',peakprom);
    % inhalation end
    [pks_in,locs_in] = findpeaks(-TH_filt,'MinPeakProminence',peakprom);
    
%     figure;
%     plot(RespirationData(:,1),TH_filt);
%     hold on
%     plot(RespirationData(locs_ex,1),pks_ex,'vk');
%     plot(RespirationData(locs_in,1),-pks_in,'.r');
    
    SniffTS = [];
    % sanity checks - missed peaks or double peaks
    for i = 1:numel(locs_in)-2
        putatives = intersect(find(locs_ex>=locs_in(i)), find(locs_ex<=locs_in(i+1)));
        if numel(putatives) > 1
            % multiple detections
            figure;
            plot(RespirationData((locs_in(i)-SampleRate):(locs_in(i+1)+SampleRate),2));
            hold on
            plot(RespirationData(locs_in(i),1),-pks_in(i),'.r');
            plot(RespirationData(locs_in(i+1),1),-pks_in(i+1),'.r');
            plot(RespirationData(locs_ex(putatives(1)),1),TH_filt(locs_ex(putatives(1)),1),'vk');
            plot(RespirationData(locs_ex(putatives(2)),1),TH_filt(locs_ex(putatives(2)),1),'vk');
            keyboard;
        elseif numel(putatives) == 0
            % missed points
            figure;
            plot(RespirationData((locs_in(i)-SampleRate):(locs_in(i+1)+SampleRate),2));
            hold on
            plot(RespirationData(locs_in(i),1),-pks_in(i),'.r');
            plot(RespirationData(locs_in(i+1),1),-pks_in(i+1),'.r');
            keyboard;
        else
            SniffTS = vertcat(SniffTS, [RespirationData(locs_ex(putatives(1)),1) RespirationData(locs_in(i+1),1) RespirationData(locs_ex(putatives(1)+1),1)]);
            % [inhalation-start inhalation-end next-inhalation]
        end
    end
    
%     % some pretty plots
%     figure;
%     [~,sortorder] = sort(SniffTS(:,3)-SniffTS(:,1));
%     SniffTS_ = SniffTS(sortorder,:);
%     line(repmat(1:size(SniffTS_,1),2,1), [(SniffTS_(:,1) - SniffTS_(:,1))'; (SniffTS_(:,3) - SniffTS_(:,1))'], 'color' ,'k');
%     hold on
%     line(repmat(1:size(SniffTS_,1),2,1), [(SniffTS_(:,1) - SniffTS_(:,1))'; (SniffTS_(:,2) - SniffTS_(:,1))'], 'color' ,'r');
%     
%     figure;
%     SniffTS_ = SniffTS;
%     line(repmat(1:size(SniffTS_,1),2,1), [(SniffTS_(:,1) - SniffTS_(:,1))'; (SniffTS_(:,3) - SniffTS_(:,1))'], 'color' ,'k');
%     hold on
%     line(repmat(1:size(SniffTS_,1),2,1), [(SniffTS_(:,1) - SniffTS_(:,1))'; (SniffTS_(:,2) - SniffTS_(:,1))'], 'color' ,'r');
%     
%     figure; 
%     figure, scatter(SniffTS(:,3) - SniffTS(:,1),SniffTS(:,2) - SniffTS(:,1), 'ok');
%     set(gca,'YLim', [0 max(get(gca,'XLim'))]);

else
    SniffTS = [];
end

end