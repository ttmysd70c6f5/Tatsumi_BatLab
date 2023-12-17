lfp_dir = 'D:\Dataset\FlyingPixels\29968\230609\ephys\20230609_124348_merged.LFP';
lfp_list = dir(fullfile(lfp_dir,'*LFP_nt*.dat')); % list of LFP .dat files
lfp_data = readTrodesExtractedDataFile(fullfile(lfp_list(1).folder, lfp_list(1).name)); % Load LFP data
nChan = length(lfp_list); % number of channels
T = length(lfp_data.fields.data); % length of recording
clear lfp_data

Fs = 1500;
Fs_sp = 30000;

%%
chan = 100;
lfpData = readTrodesExtractedDataFile(fullfile(lfp_list(chan).folder, lfp_list(chan).name)); % Load LFP data

%% original script
figure
lfp = double(lfpData.fields.data);

ax1 = subplot(3,1,1);
plot(lfp(1:Fs*60))
lfp1 = lfp - lfpmed1;
ax2 = subplot(3,1,2);
plot(lfp1(1:Fs*60))
lfp2 = lfp - median(lfp);
ax3 = subplot(3,1,3);
plot(lfp2(1:Fs*60))

linkaxes([ax1 ax2 ax3],'x')


%% 