path = '/Users/Priyanka/Desktop/LABWORK_II/Data/VideoAlignedBehaviorData/Q9/20221212_4';
for i = 35:44
    filename = fullfile(path,['Q9_camB',num2str(i-1),'.avi']);
    v = VideoReader(filename);
end

currAxes = axes;
while hasFrame(v)
    vidFrame = readFrame(v);
    image(vidFrame, 'Parent', currAxes);
    currAxes.Visible = 'off';
    pause(1/v.FrameRate);
end