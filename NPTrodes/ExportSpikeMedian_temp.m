%% Load data
spFile = dir(fullfile(spDir,'*spikeband_nt*.dat')); % list of LFP .dat files
nChan = length(spFile); % number of channels
Fs = 30000;
ext_len = Fs*5;
ts_start = Fs*10;

batch_sz = 500000; % size of each batch of devided recording

spData = readTrodesExtractedDataFile(fullfile(spFile(1).folder, spFile(1).name)); % Load LFP data
rec_len = length(spData.fields.data);
if rec_len < ts_start+ext_len
    error('The specified time window exceeds the recording length.')
end

%% Compute median for each batch
disp('Computing median...')

% rec_len = varargin{1}*varargin{2};
numBatch = ceil(ext_len/batch_sz);

spMed1 = zeros(ext_len,1);

tic
for batch = 1:numBatch-1
    fprintf('Computing median for batch %d/%d...\n',batch,numBatch)
    sp = zeros(nChan,batch_sz);
    parfor ch=1:nChan
        spData = readTrodesExtractedDataFile(fullfile(spFile(ch).folder, spFile(ch).name)); % Load LFP data
        sp(ch,:) = double(spData.fields.data(ts_start+(batch-1)*batch_sz+1:ts_start+batch*batch_sz));
    end
    spMed1((batch-1)*batch_sz+1:batch*batch_sz) = median(sp,1);
    toc
end

% final batch
batch_sz_last = ext_len - (batch-1)*batch_sz; % length of final batch

fprintf('Computing median for batch %d/%d...\n',numBatch,numBatch)
sp = zeros(nChan,batch_sz_last);

parfor ch=1:nChan

    spData = readTrodesExtractedDataFile(fullfile(spFile(ch).folder, spFile(ch).name)); % Load LFP data
    sp(ch,:) = double(spData.fields.data(ts_start+(batch-1)*batch_sz+1:ts_start+(batch-1)*batch_sz+batch_sz_last));
end
spMed1((batch-1)*batch_sz+1:end) = median(sp,1);

toc


%% Save
save(fullfile(spDir,'spmedian.mat'),'spMed1')
