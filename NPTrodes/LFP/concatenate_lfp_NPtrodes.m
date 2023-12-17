function concatenate_lfp_NPtrodes(lfp_dir)
% Make a dat file containing all channels
% Directory of lfp data: lfp_dir = 'D:\Dataset\FlyingPixels\29968\230609\ephys\20230609_124348_merged.LFP';

% load metadata
lfp_list = dir(fullfile(lfp_dir,'*LFP_nt*.dat')); % list of LFP .dat files
lfp_data = readTrodesExtractedDataFile(fullfile(lfp_list(1).folder, lfp_list(1).name)); % Load LFP data
nChan = length(lfp_list); % number of channels
T = length(lfp_data.fields.data); % length of recording
clear lfp_data

lfp = zeros(nChan,T);

for ch=1:nChan
    lfp_data = readTrodesExtractedDataFile(fullfile(lfp_list(ch).folder, lfp_list(ch).name)); % Load LFP data
    lfp(ch,:) = double(lfp_data.fields.data);
end

% Save concatenated lfp data
save(fullfile(lfp_dir,'lfp.mat'),'lfp', '-v7.3')