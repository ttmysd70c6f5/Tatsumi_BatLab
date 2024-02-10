function [] = exportAllTrodesRecordingsAllBats(data_path, trodes_path, varargin)
% Run exportTrodesRecording for recordings for all animals
% data_path is the root directory for all the data.

recs = dir(fullfile(data_path,'ephys','*.rec','*','*merged.rec')); % list of merged recording
for i = 1:length(recs)
    path_to_rec_dir = recs(i).folder;
    exportTrodesRecording(path_to_rec_dir, trodes_path, varargin{:});
end

end