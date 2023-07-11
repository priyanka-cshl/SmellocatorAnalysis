function [Frames] = GetWidefieldTimeStamps(myimagingdir, Behavior_frames)
Frames = [];

% 1: Load Frame Timestamps from the widefield folder
S = load(fullfile(myimagingdir,'params.mat'));
nBlue       = S.TS.settings(4);
nViolet     = S.TS.settings(5);

% check whether there was a gap in acquiring frames
S.TS.Blue(:,3)      = 1+ 0*S.TS.Blue(:,1);
S.TS.Violet(:,3)    = 0+ 0*S.TS.Violet(:,1);
allFrames           = vertcat(S.TS.Blue(1:nBlue,:), S.TS.Violet(1:nViolet,:));
allFrames           = sortrows(allFrames,2);
median_fs   = median(diff(allFrames(:,2)));
gaps        = find(diff(allFrames(:,2))>=1.5*median_fs);

% sanity check 1 - does the frame identity of frame 1 match on either side
if allFrames(1,3) == Behavior_frames.Behavior(1,2)
    if ~isfield(Behavior_frames,'Passive') % passive tuning was done
        gaps = size(allFrames,1);
    else
        if (numel(gaps)==1) % there is also a gap in frame timestamps
            % sanity check 2 - does the frame identity of the 1st passive frame match on either side
            if allFrames(gaps+1,3) == Behavior_frames.Passive(1,2)
                
                % get timestamps of all blue frames - converted to passive tuning timebase
                b_blue = Behavior_frames.Passive(find(Behavior_frames.Passive(:,2)),1);
                % pad any missing frames on the behavior side with NaNs
                
                w_blue = allFrames(gaps + find(allFrames(gaps+1:end,3)),2);
                
                % if there are no frame drops on the widefield side & behavior side
                if ~any(diff(b_blue)>=1.5*median(w_blue)) && ~any(diff(w_blue)>=1.5*median(w_blue))
                    % just copy over behavior timestamps
                    Frames.Passive        = Behavior_frames.Passive(find(Behavior_frames.Passive(:,2)),:);
                    Frames.Passive(:,2)   = numel(find(allFrames(1:gaps,3))) + (1:numel(b_blue)); % just frame numbers
                end
            else
                keyboard;
            end
        else
            keyboard;
        end
    end
    % get timestamps of all blue frames - converted to Behavior timebase
    b_blue = Behavior_frames.Behavior(find(Behavior_frames.Behavior(:,2)),1);
    % pad any missing frames on the behavior side with NaNs
    
    w_blue = allFrames(find(allFrames(1:gaps,3)),2);
    
    % if there are no frame drops on the widefield side & behavior side
    if ~any(diff(b_blue)>=1.5*median(w_blue)) && ~any(diff(w_blue)>=1.5*median(w_blue))
        % just copy over behavior timestamps
        Frames.Behavior        = Behavior_frames.Behavior(find(Behavior_frames.Behavior(:,2)),:);
        Frames.Behavior(:,2)   = 1:numel(b_blue); % just frame numbers
    end
    
end

end
