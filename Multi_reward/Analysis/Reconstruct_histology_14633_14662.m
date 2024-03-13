%% Reconstructing probe trajectory
% Animal: 14633, 14662
%% 14633: MEC
% Directory
ephysDir = 'Z:\users\Tatsumi\data\MultiReward\multibat\240207\ephys\20240207_103632.rec\14633\Processed';
fig_dir = 'Z:\users\Tatsumi\data\MultiReward\multibat\240207\analysis\figure\ephys';
mkdir(fig_dir)

% general info
batid = '14633';
probeNum = 1; % probe number

% Load data
lfpData = loadTrodesLFP(ephysDir);
lfp = double(lfpData.lfp{probeNum}); % LFP traces (channels x recording length)
channel_depth = lfpData.channelPositions{probeNum}(:,2); % depth of channels

%% LFP correlation
% The pearson correlation for the randomly sampled LFP traces between any channel pair is computed.
PlotLFPcorrelation(lfp,channel_depth,'spsize',30000)
% save figure
saveas(gcf,fullfile(fig_dir,sprintf('LFP_correlation_%s.png',batid)))

%% LFP power spectrogram
bandpass_freq = [2,10];
% PlotLFPspectrogram(lfp,channel_depth,bandpass_freq)
PlotLFPspectrogram(lfp,channel_depth,bandpass_freq)
% save figure
saveas(gcf,fullfile(fig_dir,sprintf('LFP_spectrogram_%s.png',batid)))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 14633: HPC
% Directory
ephysDir = 'Z:\users\Tatsumi\data\MultiReward\multibat\240206\ephys\20240206_103805.rec\14633\Processed';
fig_dir = 'Z:\users\Tatsumi\data\MultiReward\multibat\240206\analysis\figure\ephys';
mkdir(fig_dir)

% general info
batid = '14633';
probeNum = 2; % probe number

% Load data
lfpData = loadTrodesLFP(ephysDir);
lfp = double(lfpData.lfp{probeNum}); % LFP traces (channels x recording length)
channel_depth = lfpData.channelPositions{probeNum}(:,2); % depth of channels

%% LFP correlation

% The pearson correlation for the randomly sampled LFP traces between any channel pair is computed.
PlotLFPcorrelation(lfp,channel_depth,'spsize',30000)
% save figure
saveas(gcf,fullfile(fig_dir,sprintf('LFP_correlation_%s_probe%d.png',batid,probeNum)))

%% LFP power spectrogram
bandpass_freq = [2,10];
PlotLFPspectrogram(lfp,channel_depth,bandpass_freq)
% save figure
saveas(gcf,fullfile(fig_dir,sprintf('LFP_spectrogram_%s_probe%d.png',batid,probeNum)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 14662: HPC-MEC
% Directory
ephysDir = 'Z:\users\Tatsumi\data\MultiReward\multibat\240209_1\ephys\20240209_105848.rec\14662\Processed';
fig_dir = 'Z:\users\Tatsumi\data\MultiReward\multibat\240209_1\analysis\figure\ephys';
mkdir(fig_dir)

% general info
batid = '14662';
probeNum = 1; % probe number

% Load data
lfpData = loadTrodesLFP(ephysDir);
lfp = double(lfpData.lfp{probeNum}); % LFP traces (channels x recording length)
channel_depth = lfpData.channelPositions{probeNum}(:,2); % depth of channels

%% LFP correlation

% The pearson correlation for the randomly sampled LFP traces between any channel pair is computed.
PlotLFPcorrelation(lfp,channel_depth,'spsize',30000)
% save figure
saveas(gcf,fullfile(fig_dir,sprintf('LFP_correlation_%s_probe%d.png',batid,probeNum)))

%% LFP power spectrogram
bandpass_freq = [2,10];
PlotLFPspectrogram(lfp,channel_depth,bandpass_freq)
% save figure
saveas(gcf,fullfile(fig_dir,sprintf('LFP_spectrogram_%s_probe%d.png',batid,probeNum)))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 14662: HPC
% Directory
ephysDir = 'Z:\users\Tatsumi\data\MultiReward\multibat\240209_2\ephys\20240209_154525.rec\14662\Processed';
fig_dir = 'Z:\users\Tatsumi\data\MultiReward\multibat\240209_2\analysis\figure\ephys';
mkdir(fig_dir)

% general info
batid = '14662';
probeNum = 1; % probe number

% Load data
lfpData = loadTrodesLFP(ephysDir);
lfp = double(lfpData.lfp{probeNum}); % LFP traces (channels x recording length)
channel_depth = lfpData.channelPositions{probeNum}(:,2); % depth of channels

%% LFP correlation

% The pearson correlation for the randomly sampled LFP traces between any channel pair is computed.
PlotLFPcorrelation(lfp,channel_depth,'spsize',30000)
% save figure
saveas(gcf,fullfile(fig_dir,sprintf('LFP_correlation_%s_probe%d.png',batid,probeNum)))

%% LFP power spectrogram
bandpass_freq = [2,10];
PlotLFPspectrogram(lfp,channel_depth,bandpass_freq)
% save figure
saveas(gcf,fullfile(fig_dir,sprintf('LFP_spectrogram_%s_probe%d.png',batid,probeNum)))

%% Putative ripples
fig = figure;
fontsize(fig, 24, "points")
set(gcf, 'units','normalized','Position',[0 0 1 1]);
tiledlayout(5,10,"TileSpacing",'compact')
for i = 1:50
    nexttile
    ExtractRipples(lfp(i*7,1:600*1500),true);
end