function [MyData, DataTags] = OdorLocationSanityCheck(MyData, DataTags)

global errorflags; % [digital-analog sample drops, timestamp drops, RE voltage drift, motor slips]

checkflags = [0 0 0];

MotorCol = find(cellfun(@isempty,regexp(DataTags,'Motor'))==0);
EncoderCol = find(cellfun(@isempty,regexp(DataTags,'Encoder'))==0);

% check that TEENSY analog output had no drift
[coeff,gof] = fit(MyData(:,MotorCol),MyData(:,EncoderCol),'poly1');
checkflags(1) = (gof.rsquare > 0.98);

% check if Motor Position values were ~0 whenever the Motor physically
% crossed home - interrupted the photodiode
HomeCol = find(cellfun(@isempty,regexp(DataTags,'HomeSensor'))==0);
foo = MyData(:,MotorCol); 
foo(find(MyData(:,HomeCol)==0)) = NaN;
fooLims = median(foo,'omitnan') + [std(foo,'omitnan') -std(foo,'omitnan')];
checkflags(2) = ~any(abs(fooLims)>4);

% figure;
% scatter(MyData(:,MotorCol),MyData(:,EncoderCol))
% hold on
% scatter(foo,MyData(:,EncoderCol),'r')

% check if Motor Position values were +/- 8 whenever the lever was supposed
% to be in the TargetZone
TZCol = find(cellfun(@isempty,regexp(DataTags,'InTargetZone'))==0);
foo = MyData(:,MotorCol); 
foo(find(MyData(:,TZCol)==0)) = NaN;
fooLims = median(foo,'omitnan') + [std(foo,'omitnan') -std(foo,'omitnan')];
checkflags(3) = ~any(abs(fooLims)>10);

% scatter(foo,MyData(:,EncoderCol),'.k')

if ~any(checkflags)
     disp('Warning: session failed sanity check');
else
    disp('session passed sanity check');
    % delete the Encoder Col and Homesensor column
    MyData(:,[EncoderCol HomeCol]) = [];
    DataTags([EncoderCol HomeCol],:) = [];
end

errorflags(3:5) = checkflags;
end