function newEphysDir = GetEphysDir(ephys_rootdir,batid)
% function to get the path of the recording directory for the specified bat

newEphysDir = fullfile(ephys_rootdir,batid,'Processed');

end