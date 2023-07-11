function [Frame_TS] = GetPCOFrameTriggers(MyData)

%% Get Camera Timestamps for the closed loop stretch
C_on        = MyData((find(diff(MyData(:,17))==-1) + 1),1); % blue frames
C_on(:,2)   = 1+ 0*C_on(:,1);
C_off       = MyData((find(diff(MyData(:,17))==+1) + 1),1); % violet frames
C_off(:,2)  = 0+ 0*C_off(:,1);

Frame_TS = vertcat(C_off,C_on);
Frame_TS = sortrows(Frame_TS, 1);
median_fs = median(diff(Frame_TS(:,1)));

frame_drops = any(diff(Frame_TS(:,1))>=1.5*median_fs);
if ~frame_drops
    disp('no frames dropped during behavior');
else
    disp('warning: frames may have been dropped during behavior');
    while frame_drops
        idx = find(diff(abs(Frame_TS(:,1)))>=1.5*median_fs,1,'first') + 1;
        
        % just use the next frame to be safe - last transition might be fake
        frame_diff = Frame_TS(idx+1,1) - Frame_TS(idx-1,1);
        missing_frames = round((frame_diff - median_fs)/median_fs);
        % sanity check
        if ~mod(missing_frames,2) && isequal(Frame_TS(idx-1,2),Frame_TS(idx+1,2)) % even frames
            disp('unresolved frame drops');
            keyboard;
        elseif mod(missing_frames,2) && ~isequal(Frame_TS(idx-1,2),Frame_TS(idx+1,2)) % odd frames
            disp('unresolved frame drops');
            keyboard;
        else
            pad = mod(Frame_TS(idx-1,2)+(1:missing_frames),2);
            ts  = linspace(Frame_TS(idx-1,1),Frame_TS(idx+1,1),(missing_frames + 2));
            Frame_TS = vertcat(Frame_TS(1:idx-1,:),...
                               [-ts(2:end-1)' pad'],...
                                Frame_TS((idx+1):end,:));
        end
        
        frame_drops = any(diff(abs(Frame_TS(:,1)))>=1.5*median_fs);
    end
end

% another sanity check
if ~(numel(unique(Frame_TS(2:2:end,2))) == 1) ...
        || ~(numel(unique(Frame_TS(1:2:end,2))) == 1)
    disp('unresolved frame drops');
    keyboard;
end
