function convert_Trodes_dat2binary(data_folder,output_folder,stream,gain,convert_uv,varargin)
% stream: LFP or spikeband
% samplePerBatch: batch size
% convert_uv: True -> convert values to uv

current_folder = pwd;
cd(data_folder)

%% Convert the data in LFP band to binary file
% % Output folder
[~,foldername,~] = fileparts(pwd);
% output_folder = fullfile(folderpath,[foldername,'.binary']);
% if ~isfolder(output_folder)
%     mkdir(output_folder)
% end

% Convert to uV
% gain = 50;
if convert_uv
    scale_to_uv = (600/32767)*(1000/gain);
else
    scale_to_uv = 1;
end

% Number of probes
num_chans = zeros(1,3);
for probe = 1:3
    num_chans(probe) = length(dir(fullfile(sprintf('*%s_nt%d*',stream,probe))));
end
num_probe = find(num_chans~=0,1,"last");

% Make a list of all dat files
stream_file = struct([]);
for probe = 1:num_probe
    stream_file(probe).filepath = dir(sprintf('*%s_nt%d*',stream,probe));
end

% timestamps
ts_file = dir('*timestamps.dat');
ts_data = readTrodesExtractedDataFile(fullfile(ts_file.folder,ts_file.name)); % Load timestamps for LFP
ts_lfp = ts_data.fields.data; % time (point)

%% Convert Spikegadgets .dat to non-header binary data
if nargin>5
    samplePerBatch = varargin{:};
    num_batch = ceil(length(ts_lfp)/samplePerBatch);
    batch_size = zeros(1,num_batch);
    for batch = 1:num_batch
        if batch < num_batch
            batch_size(batch) = samplePerBatch;
        else
            batch_size(batch) = rem(length(ts_lfp),samplePerBatch);
        end
    end
else
    samplePerBatch = length(ts_lfp);
    num_batch = ceil(length(ts_lfp)/samplePerBatch);
    batch_size = samplePerBatch;
end

for probe = 1:num_probe
    offset = 0;
    fid = fopen(fullfile(output_folder,sprintf('%s.%s_probe%d.dat',foldername,stream,probe)),'w');
    for batch = 1:num_batch
        offset = offset+(batch-1)*samplePerBatch;
        stream_batch = int16(zeros(num_chans(probe),batch_size(batch)));
        for ch = 1:num_chans(probe)
            current_file = stream_file(probe).filepath(ch);
            stream_data = readTrodesExtractedDataFile_new(fullfile(current_file.folder,current_file.name),offset,batch_size(batch));
            stream_batch(ch,:) = stream_data.fields.data*scale_to_uv;
        end
        % write to a binary file
        stream_batch = reshape(stream_batch,numel(stream_batch),1);
        fwrite(fid,stream_batch,'int16');
    end
    fclose(fid);
end

%% Read the binary file
num_chan = 200;

% Original file
ap_dir = 'C:\Users\Tatsumi\Documents\Data\KQTY_NP\32623\20230907_1min\20230907_161730_merged_split1.LFP';
ap_file = dir(fullfile(ap_dir,'*nt2*.dat'));
% samples = 100;
fname = fullfile(ap_file(1).folder,ap_file(1).name);
D = readTrodesExtractedDataFile(fname);
samples = length(D.fields.data);
D_original = zeros(num_chan,samples);
for ch = 1:length(ap_file)
    fname = fullfile(ap_file(ch).folder,ap_file(ch).name);
    D = readTrodesExtractedDataFile(fname);
    D_original(ch,:) = double(D.fields.data);
end
%%
% concatenated file
fname = 'C:\Users\Tatsumi\Documents\Data\KQTY_NP\32623\20230907_1min\20230907_161730_merged_split1.binary\20230907_161730_merged_split1.LFP_mmap_probe2.dat';
binary_file = dir(fname);
samples = binary_file.bytes/2/num_chan;
fid = fopen(fname,'rb');
D_new = fread(fid,[num_chan,samples],'int16');
fclose(fid);
% % Original file
% fname = 'C:\Users\Tatsumi\Documents\Data\KQTY_NP\32623\20230907_161730.rec\20230907_161730_merged.spikeband';
% binary_file = dir(fname);
% samples = binary_file.bytes/2/num_chan;
% D2.Data = memmapfile(fname,'Format',{'int16' [num_chan samples] 'mapped'});
% D2 = D2.Data.Data.mapped;
%%
cd(current_folder)