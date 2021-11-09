function rgb= Plot_Colors(string)

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
end


