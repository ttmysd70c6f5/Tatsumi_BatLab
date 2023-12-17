%% Examine the contents of kilosort .dat file
%% Load kilosort file
filename = dir(fullfile(pwd,'*kilosort','*merged.probe1.dat'));
data_file = fullfile(filename(1).folder,filename(1).name);
fid = fopen(data_file,'rb');

%% Read header
nchan = 192;
sample_to_read = 100000;
head_values = fread(fid,[nchan,sample_to_read],"int16");
frewind(fid)

amp_chan = max(abs(head_values),[],2);
figure; plot(amp_chan)
%% Close file
fclose(fid)


%% Compare non-filtered spikeband and kilosort .dat file
nchan = 192;
sample_to_read = 100000;
% spikeband
data_file = dir(fullfile(pwd,'*spikeband','*nt1*'));
head_ap = zeros(nchan,sample_to_read);
for ch = 1:length(data_file)
    temp_data = readTrodesExtractedDataFile(fullfile(data_file(ch).folder,data_file(ch).name));
    head_ap(ch,:) = temp_data.fields.data(1:sample_to_read);
end
% kilosort
filename = dir(fullfile(pwd,'*kilosort','*probe1.dat'));
data_file = fullfile(filename(1).folder,filename(1).name);
fid = fopen(data_file,'rb');
head_ks = fread(fid,[nchan,sample_to_read],"int16");
fclose(fid);

%% Examine the contents of raw .dat file
headerText = fread(fid,1000,"uint8");
headerText = char(headerText)

%% Find the fist value
fieldsLoc  = strfind(headerText,'First_timestamp');
fseek(fid, fieldsLoc+25, -1);

% Read the first 1500 columns
values = fread(fid,[200,1500],'uint16=>int16');

%%
fname = dir(fullfile(pwd,'*LFP','*nt1001*'));
fname = fullfile(fname.folder,fname.name)
lfp_data = readTrodesExtractedDataFile(fname);

%% Check header format of spikegadgets data
% file_path = 'C:\Users\Tatsumi\Documents\GitHub\SpikeInterface\Tutorial\spikegadgets_tutorial.rec';
file_path = 'C:\Users\Tatsumi\Documents\Data\KQTY_NP\32623\20230907_161730.rec\20230907_161730_merged_split1.rec';
fid = fopen(file_path,'rb');
header_1 = fread(fid,100000,"uint8");
header_1 = char(header_1');
frewind(fid)


%%
header_1 = char(fread(fid,500000,"uint8")');
% numbytes_loc = strfind(header_1,'numBytes')
device_loc = strfind(header_1,'name')

for i = 1:length(device_loc)
    % fseek(fid,numbytes_loc(i)-100,-1);
    fseek(fid,device_loc(i)-100,-1);
    char(fread(fid,200,"uint8")')
end
frewind(fid)