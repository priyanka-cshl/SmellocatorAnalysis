function [Frame_TS] = GetPCOFrameTriggers(MyData)

%% Get Camera Timestamps for the closed loop stretch
C_on        = MyData((find(diff(MyData(:,17))==-1) + 1),1); % blue frames
C_on(:,2)   = 1+ 0*C_on(:,1);
C_off       = MyData((find(diff(MyData(:,17))==+1) + 1),1); % violet frames
C_off(:,2)  = 0+ 0*C_off(:,1);

Frame_TS = vertcat(C_off,C_on);
Frame_TS = sortrows(Frame_TS, 1);
median_fs = median(diff(Frame_TS(:,1)));

if ~any(diff(Frame_TS(:,1))>=1.5*median_fs)
    disp('no frames dropped during behavior');
else
    MyFrames = Frame_TS(1,:);
    disp('warning: frames may have been dropped during behavior');
    i = 2;
    while i<=size(Frame_TS,1)
        frame_diff = Frame_TS(i,1) - Frame_TS(i-1,1);
        if frame_diff >= 1.5*median_fs
            % just use the next frame to be safe - last transition might be
            % fake
            frame_diff = Frame_TS(i+1,1) - Frame_TS(i-1,1);
            missing_frames = round((frame_diff - median_fs)/median_fs);
            % sanity check
            if ~mod(missing_frames,2) && isequal(Frame_TS(i-1,2),Frame_TS(i+1,2)) % even frames
                disp('unresolved frame drops');
                keyboard;
            elseif mod(missing_frames,2) && ~isequal(Frame_TS(i-1,2),Frame_TS(i+1,2)) % odd frames
                disp('unresolved frame drops');
                keyboard;
            else
                pad = mod(MyFrames(end,2)+(1:missing_frames),2);
                MyFrames = vertcat(MyFrames,[NaN*pad' pad'],Frame_TS(i+1,:));
                i = i + 2;
            end
        else
            MyFrames = vertcat(MyFrames,Frame_TS(i,:));
            i = i + 1;
        end
    end
    Frame_TS = MyFrames;
end

% another sanity check
if ~(numel(unique(Frame_TS(2:2:end,2))) == 1) ...
        || ~(numel(unique(Frame_TS(1:2:end,2))) == 1)
    disp('unresolved frame drops');
    keyboard;
end
