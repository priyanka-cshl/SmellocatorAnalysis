function [Paths] = WhichComputer()

% read computer name
[~,computername] = system('hostname');
computername = deblank(computername);

% default paths
Paths.Code                      = '/opt'; % where all code files are
Paths.Grid.Behavior             = '/mnt/grid-hs/pgupta/Behavior'; % raw behavior, tuning mtalab files
Paths.Grid.Ephys{1}             = '/mnt/grid-hs/pgupta/EphysData'; % raw oeps files - for PCX batch
Paths.Grid.Ephys{2}             = '/mnt/grid-hs/mdussauz/ephysdata/lever_task/Batch'; % raw oeps files - for batch O,MO, J
Paths.Local.Ephys_processed     = '/mnt/data/Sorted'; % local copy where sorted, curated spike data is stored
Paths.Grid.Ephys_processed      = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Ephys'; % local copy where sorted, curated spike data is stored
Paths.Local.Behavior_processed  = '/mnt/data/Behavior'; % local copy where sorted, curated spike data is stored
Paths.Grid.Behavior_processed   = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior'; % local copy where sorted, curated spike data is stored

switch computername
    case 'andaman'
        % use defaults
end

end