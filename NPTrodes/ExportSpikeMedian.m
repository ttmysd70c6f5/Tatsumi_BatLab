function ExportSpikeMedian(spDir,varargin)

%% Load data
spFile = dir(fullfile(spDir,'*spikeband_nt*.dat')); % list of LFP .dat files
nChan = length(spFile); % number of channels


batch_sz = 100000; % size of each batch of devided recording

%% Compute median for each batch
disp('Computing median...')

if nargin == 1
    spData = readTrodesExtractedDataFile(fullfile(spFile(1).folder, spFile(1).name)); % Load LFP data
    rec_len = length(spData.fields.data); % length of recording
    numBatch = ceil(rec_len/batch_sz);
    clear spData
    

    spMed1 = zeros(rec_len,1);
    for batch = 1:numBatch-1
        fprintf('Computing median for batch %d/%d...\n',batch,numBatch)
        sp = zeros(nChan,batch_sz);
        parfor ch=1:nChan
            spData = readTrodesExtractedDataFile(fullfile(spFile(ch).folder, spFile(ch).name)); % Load LFP data
            sp(ch,:) = double(spData.fields.data((batch-1)*batch_sz+1:batch*batch_sz));
        end
        spMed1((batch-1)*batch_sz+1:batch*batch_sz) = median(sp,1);
    end
    
    % final batch
    batch_sz_last = rec_len - (batch-1)*batch_sz; % length of final batch
    
    fprintf('Computing median for batch %d/%d...\n',numBatch,numBatch)
    sp = zeros(nChan,batch_sz_last);
    for ch=1:nChan
        spData = readTrodesExtractedDataFile(fullfile(spFile(ch).folder, spFile(ch).name)); % Load LFP data
        sp(ch,:) = double(spData.fields.data((batch-1)*batch_sz+1:end));
    end
    spMed1((batch-1)*batch_sz+1:end) = median(sp,1);

elseif nargin == 2
    rec_len = varargin{1}*varargin{2};
    numBatch = ceil(rec_len/batch_sz);

    spMed1 = zeros(rec_len,1);
    
    for batch = 1:numBatch-1
        fprintf('Computing median for batch %d/%d...\n',batch,numBatch)
        sp = zeros(nChan,batch_sz);
        cnt = 0;
        parfor ch=1:nChan
            cht = cnt + 1;
            fprintf('Processed %d/%d electrodes...\n',cnt,nChan)
            spData = readTrodesExtractedDataFile(fullfile(spFile(ch).folder, spFile(ch).name)); % Load LFP data
            sp(ch,:) = double(spData.fields.data((batch-1)*batch_sz+1:batch*batch_sz));
        end
        spMed1((batch-1)*batch_sz+1:batch*batch_sz) = median(sp,1);
    end
    
    % final batch
    batch_sz_last = rec_len - (batch-1)*batch_sz; % length of final batch
    
    fprintf('Computing median for batch %d/%d...\n',numBatch,numBatch)
    sp = zeros(nChan,batch_sz_last);
    cnt = 0;
    parfor ch=1:nChan
        cht = cnt + 1;
        fprintf('Processed %d/%d electrodes...\n',cnt,nChan)
        spData = readTrodesExtractedDataFile(fullfile(spFile(ch).folder, spFile(ch).name)); % Load LFP data
        sp(ch,:) = double(spData.fields.data((batch-1)*batch_sz+1:(batch-1)*batch_sz+batch_sz_last));
    end
    spMed1((batch-1)*batch_sz+1:end) = median(sp,1);

end

%% Save
save(fullfile(spDir,'spmedian.mat'),'spMed1')
