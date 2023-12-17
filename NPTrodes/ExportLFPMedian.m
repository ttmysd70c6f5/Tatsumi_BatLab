function ExportLFPMedian(lfp_dir)

%% Load data
lfp_file = dir(fullfile(lfp_dir,'*LFP_nt*.dat')); % list of LFP .dat files
nChan = length(lfp_file); % number of channels

lfp_data = readTrodesExtractedDataFile(fullfile(lfp_file(1).folder, lfp_file(1).name)); % Load LFP data
rec_len = length(lfp_data.fields.data); % length of recording
clear lfp_data

batch_sz = 1000000; % size of each batch of devided recording
numBatch = ceil(rec_len/batch_sz);

%% Compute median for each batch
disp('Computing median...')
lfpmed1 = zeros(rec_len,1);
for batch = 1:numBatch-1
    fprintf('Computing median for batch %d/%d...\n',batch,numBatch)
    lfp = zeros(nChan,batch_sz);
    for ch=1:nChan
        lfp_data = readTrodesExtractedDataFile(fullfile(lfp_file(ch).folder, lfp_file(ch).name)); % Load LFP data
        lfp(ch,:) = double(lfp_data.fields.data((batch-1)*batch_sz+1:batch*batch_sz));
    end
    lfpmed1((batch-1)*batch_sz+1:batch*batch_sz) = median(lfp,1);
end

% final batch
batch_sz_last = rec_len - (batch-1)*batch_sz; % length of final batch

fprintf('Computing median for batch %d/%d...\n',numBatch,numBatch)
lfp = zeros(nChan,batch_sz_last);
for ch=1:nChan
    lfp_data = readTrodesExtractedDataFile(fullfile(lfp_file(ch).folder, lfp_file(ch).name)); % Load LFP data
    lfp(ch,:) = double(lfp_data.fields.data((batch-1)*batch_sz+1:end));
end
lfpmed1((batch-1)*batch_sz+1:end) = median(lfp,1);

%% Save
save(fullfile(lfp_dir,'lfpmedian.mat'),'lfpmed1')
