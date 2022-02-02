MySession = '/mnt/data/Processed/Behavior/O3/O3_20211005_r0_processed.mat'; % session path - leave empty to get browser pop up
LoadProcessedSession; % loads relevant variables

% create a corresponding timestamp vector
Timestamps = TracesOut(:,find(strcmp(ColNames,'Timestamps')))'; % in behavior timebase

N = size(SingleUnits,2);

%% sort units by tetrode - to match session viewer
clear foo
for i = 1:N
    foo(i,:) = [SingleUnits(i).tetrode SingleUnits(i).id];
end
[~,SortedByTetrodes] = sort(foo(:,1));
UnitOrder = SortedByTetrodes;

%%
ChosenUnits = [58 35 34 55 21];
STAwindow = [-1 1];
snippetlength = diff(STAwindow)*SampleRate;
%F = fieldnames(TracesOut);
% for every unit
for i = 1:numel(ChosenUnits)
    MyUnit = ChosenUnits(i);
    % for every spike (in behavior timebase
        thisunitspiketimes = SingleUnits(MyUnit).spikes - TimestampAdjuster + STAwindow(1);
        % ignore spikes that precede behavior start
        thisunitspiketimes(thisunitspiketimes<0) = [];
        % ignore spikes that follow behavior stop
        thisunitspiketimes(thisunitspiketimes>=(Timestamps(end)-2*STAwindow(2))) = [];
        window = [];
        Snippets = NaN*zeros(snippetlength+1,numel(ColNames),numel(thisunitspiketimes));
        for j = 1:numel(thisunitspiketimes)
            idx1 = find(Timestamps>=thisunitspiketimes(j),1,'first');
            idx2 = idx1 + snippetlength;
%            for k = 1:numel(F)-1
                Snippets(:,:,j) = TracesOut(idx1:idx2,:);
%                 mysnippet = TracesOut.(F{k}){1}(idx1:idx2);
%                 if j == 1
%                     Snippets(:,k) = mysnippet;
%                 else
%                     Snippets(:,k) = Snippets(:,k) + mysnippet;
%                 end
%            end
        end
        % get the average trace
        STA(:,:,i) = mean(Snippets,3,'omitnan');
        %STASTD(:,:,i) = std(Snippets,3,'omitnan');
%         for k = 1:numel(F)-1
%             STA(:,k,i) = Snippets
%         end
end
    
% TraceNames = {'Lever' 'Motor' 'Sniffs' 'Licks' 'Rewards' 'Trial' 'TargetZone'};
%% plotting
figure
for i = 1:5
    subplot(2,5,i);
    plot(STA(:,8:10,i));
    hold on
    line( SampleRate*[1 1],[0 50],'Color','k','LineStyle',':');
    line( SampleRate*[0.5 0.5],[0 50],'Color','k','LineStyle',':');
    
    subplot(2,5,i+5);
    plot(STA(:,1,i),'k');
    hold on
    line( SampleRate*[1 1],[0 5],'Color','k','LineStyle',':');
    line( SampleRate*[0.5 0.5],[0 5],'Color','k','LineStyle',':');
end