function LFP_CAR(lfp_dir)

%% Load data
lfp_list = dir(fullfile(lfp_dir,'*LFP_nt*.dat')); % list of LFP .dat files
lfp_data = readTrodesExtractedDataFile(fullfile(lfp_list(1).folder, lfp_list(1).name)); % Load LFP data
nChan = length(lfp_list); % number of channels
rec_len = length(lfp_data.fields.data); % length of recording
clear lfp_data

batch_sz = 1000000; % size of each batch of devided recording
numBatch = floor(rec_len/batch_sz);

%% Compute median for each batch
disp('Computing median...')
lfp_car = zeros(rec_len,1);
for batch = 1:numBatch-1
    fprintf('Computing median for batch %d/%d...\n',batch,numBatch)
    lfp = zeros(nChan,batch_sz);
    for ch=1:nChan
        lfp_data = readTrodesExtractedDataFile(fullfile(lfp_list(ch).folder, lfp_list(ch).name)); % Load LFP data
        lfp(ch,:) = double(lfp_data.fields.data((batch-1)*batch_sz+1:batch*batch_sz));
    end
    lfp = lfp - median(lfp,2);
    lfp = lfp - median(lfp,1);
    lfp_car((batch-1)*batch_sz+1:batch*batch_sz) = lfp;
end

% final batch
batch_sz_last = rec_len - (numBatch-1)*batch_sz; % length of final batch

fprintf('Computing median for batch %d/%d...\n',numBatch,numBatch)
lfp = zeros(nChan,batch_sz_last);
for ch=1:nChan
    lfp_data = readTrodesExtractedDataFile(fullfile(lfp_list(ch).folder, lfp_list(ch).name)); % Load LFP data
    lfp(ch,:) = double(lfp_data.fields.data((numBatch-1)*batch_sz+1:end));
end
lfp = lfp - median(lfp,2);
lfp = lfp - median(lfp,1);
lfp_car((numBatch-1)*batch_sz+1:end) = median(lfp,1);

%% Save
save(fullfile(lfp_dir,'lfp_CAR.mat'),'lfp_car')
