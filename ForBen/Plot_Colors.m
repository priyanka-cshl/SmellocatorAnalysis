function colorout = Plot_Colors(string, extras)

if nargin<2
    extras = [];
end

switch string
    case 'g'
        rgb=[78 154 6]/256;
    case 'b'
        rgb=[52 101 164]./256;
    case 'pd'%purple light
        rgb=[92 53 102]./256;
    case 'pl'%purple dark
        rgb=[173 127 168]./256;
    case 'o'%orange
        rgb=[245 121 0]./256;
    case 'r'%red
        rgb=[239 41 41]./256;
    case 't'%turquize
        rgb=[0 139 139]./256;
    case 'p'%pink
        rgb=[199 21 133]./256;
    case 'k'
        rgb = [0 0 0];
    case 'Odor1'
        rgb = [.8 .8 .8];
    case 'Odor2'
        rgb = [0.8941    0.9412    0.9020];
    case 'Odor3'
        rgb = [0.8706    0.9216    0.9804];
    case 'Odor4'
        rgb = [0.93    0.84    0.84];
    case 'TZ'
        rgb = [1 1 0];
    case 'Paletton'
        % blues
        rgb(1,:,1) = [5.9 25.9 30.6]/100;
        rgb(1,:,2) = [14.5 36.1 41.2]/100;
        rgb(1,:,3) = [26.7 47.1 52.2]/100;
        rgb(1,:,4) = [42.7 59.6 63.5]/100;
        % reds
        rgb(2,:,1) = [49.4 12.9 8.2]/100;
        rgb(2,:,2) = [66.7 27.5 22.4]/100;
        rgb(2,:,3) = [83.9 46.7 42]/100;
        rgb(2,:,4) = [100 70.2 66.3]/100;
        % yellows
        rgb(3,:,1) = [49.4 36.9 8.2]/100;
        rgb(3,:,2) = [66.7 52.9 22.4]/100;
        rgb(3,:,3) = [83.9 71.74 42]/100;
        rgb(3,:,4) = [100 89.8 66.3]/100;
   
        
end

if ~isempty(extras)
    colorout = squeeze(rgb(extras(1),:,extras(2)));
else
    colorout = rgb;
end


