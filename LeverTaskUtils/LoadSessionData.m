function [MyData, MyParams, DataTags, varout1, varout2] = LoadSessionData(FileName, TuningFlag, PIDflag)

%% load the file
Temp = load(FileName,'session_data');
MyTraces = Temp.session_data.trace; % data
MyParams = Temp.session_data.params; % settings
[nrows, ~] = size(MyTraces);
if ~TuningFlag
    varout1 = Temp.session_data.ForNextSession(8); % 0 = from left, 1 = from right
    varout2 = [];
else
    varout1 = Temp.session_data.TrialSequence; % location, odor
    if isfield(Temp.session_data,'LocationSequence')
        varout2 = Temp.session_data.LocationSequence;
        if (size(varout2,1)==(size(varout1,1)-1)) && (varout1(end,1)>997)
            % if last trial happens to be a passive replay then there is no
            % entry in the location sequence for that - it should just be
            % all zeros
            varout2(end+1,:) = 0;
        end
    else
        varout2 = zeros(size(varout1,1),1);
    end
end
DataTags = {'Lever'; ...
            'Encoder'; ...
            'TrialON'; ...
            'InTargetZone'; ...
            'InRewardZone'; ...
            'Rewards'; ...
            'Licks'};

if PIDflag == 1
    DataTags{1} = 'PID';
    MyData = MyTraces(:,[6 3 7:11]);
else
    MyData = MyTraces(:,[1 3 7:11]);
    
    % add 3 cols in the beginning [timestamps TZoneHighLim TZoneLowLim ...
    % add 2 cols in the end [WhichPerturbation PerturbationValue]
    MyData = horzcat(Temp.session_data.timestamps, zeros(nrows,2), MyData, zeros(nrows,2));
    DataTags = cat(1,{'Timestamps'; 'TZoneHighLim'; 'TZoneLowLim'}, DataTags(:), {'WhichPerturbation'; 'PerturbationValue'});
    
    % append motor location
    MyData(:,13) = MyTraces(:,4);
    DataTags = cat(1, DataTags(:), {'Motor'});
    
    if find(ismember(Temp.session_data.trace_legend,'homesensor'))
        whichcol = find(ismember(Temp.session_data.trace_legend,'homesensor'));
        MyData(:,14) = MyTraces(:,whichcol);
        DataTags = cat(1, DataTags(:), {'HomeSensor'});
    end
    
    if find(ismember(Temp.session_data.trace_legend,'respiration'))
        whichcol = find(ismember(Temp.session_data.trace_legend,'respiration'));
        MyData(:,15) = MyTraces(:,whichcol);
        DataTags = cat(1, DataTags(:), {'respiration'});
    end
    
    if find(ismember(Temp.session_data.trace_legend,'thermistor'))
        whichcol = find(ismember(Temp.session_data.trace_legend,'thermistor'));
        MyData(:,15) = MyTraces(:,whichcol);
        DataTags = cat(1, DataTags(:), {'thermistor'});
    end
    
    if find(ismember(Temp.session_data.trace_legend,'camerasync'))
        whichcol = find(ismember(Temp.session_data.trace_legend,'camerasync'));
        MyData(:,16) = MyTraces(:,whichcol);
        DataTags = cat(1, DataTags(:), {'Pgrey1'});
    end
    
    if size(MyTraces,2)>whichcol % if there was a second camera
        whichcol = whichcol + 1;
        MyData(:,17) = MyTraces(:,whichcol);
        DataTags = cat(1, DataTags(:), {'Pgrey2'});
    end
    
    if PIDflag == 2
        % replace column 10 - licks - with lick piezo
        if find(ismember(Temp.session_data.trace_legend,'lickpiezo'))
            whichcol = find(ismember(Temp.session_data.trace_legend,'lickpiezo'));
            MyData(:,10) = MyTraces(:,whichcol);
            %DataTags = cat(1, DataTags(:), {'Piezo'});
        end
    end
    
end