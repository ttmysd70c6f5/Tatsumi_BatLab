ephysDir_2probe = 'Z:\users\Tatsumi\data\MultiReward\tethered\14662\240208\ephys\20240208_155244.rec\Processed';
lfpData = loadTrodesLFP(ephysDir_2probe);

%%
% HPC-MEC
probeNum = 1;
lfp_probe1 = double(lfpData.lfp{probeNum}); % LFP traces (channels x recording length)
channel_depth = lfpData.channelPositions{probeNum}(:,2); % depth of channels

lfp_probe1 = lfp_probe1 - median(lfp_probe1,2,'omitnan');
lfp_probe1 = lfp_probe1 - median(lfp_probe1,1,'omitnan');

% sort data
[~, sorted_idx] = sort(channel_depth); % sorted channel depth and index
sorted_lfp_probe1 = lfp_probe1(sorted_idx,:); % sorted LFP

%%
% dCA1_vCA1
probeNum = 2;
lfp_probe2 = double(lfpData.lfp{probeNum}); % LFP traces (channels x recording length)
channel_depth = lfpData.channelPositions{probeNum}(:,2); % depth of channels

lfp_probe2 = lfp_probe2 - median(lfp_probe2,2,'omitnan');
lfp_probe2 = lfp_probe2 - median(lfp_probe2,1,'omitnan');

% sort data
[~, sorted_idx] = sort(channel_depth); % sorted channel depth and index
sorted_lfp_probe2 = lfp_probe2(sorted_idx,:); % sorted LFP

%%
spsize = 30000;
nChan = size(lfp_probe1,1); % number of channels
len_rec = size(lfp_probe1,2); % length of recording
fs = 1500; % sampling rate of lfp
idx = randperm(len_rec, spsize); % sampled index

lfp_merged = [sorted_lfp_probe1; sorted_lfp_probe2];
% lfp_merged = normalize(lfp_merged,2,"zscore");
lfp_corr_merged = corr(lfp_merged(:,idx)');

%%
figure
imagesc(lfp_corr_merged); % long channelmap
set(gca,'YDir','normal', 'XDir','normal', 'TickLength', [0,0], ...
    'fontsize', 10, 'DataAspectRatio', [1 1 1], ...
    'GridAlpha', 1, 'MinorGridAlpha', 1)

colormap jet
% colorbar
title('LFP correlation')
xlabel('Depth (um)')
ylabel('Depth (um)')