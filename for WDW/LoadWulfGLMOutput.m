function [PSTHOut,TimestampsOut] = LoadWulfGLMOutput(myPath,nUnits,TimestampsIn)

load(myPath);

interpolate = 0; 
if nargin == 3
    if ~isempty(TimestampsIn)
        interpolate = 1;
    end
end

Timestamps = full_active_y_pred_times;
if ~interpolate
    TimestampsOut = Timestamps;
else
    TimestampsOut = TimestampsIn;
    scaling = mean(diff(TimestampsIn))/mean(diff(Timestamps));
end
for i = 1:nUnits
    if ~interpolate
        PSTHOut(i,:) = eval(['neuron_',num2str(i-1),'.full_active_y_pred']);
    else
        myPSTH = eval(['neuron_',num2str(i-1),'.full_active_y_pred']);
        PSTHOut(i,:) = scaling*interp1(Timestamps,myPSTH,TimestampsIn);
    end
end

end