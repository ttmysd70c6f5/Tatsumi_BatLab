function lfpCorrAllChPlot(lfp_dir)

%% Set the lfp directory
% lfp_dir = 'C:\Users\Tatsumi\Documents\KQTY_NP\FlyingPixels\29968\230531\ephys\20230531_115208_merged.LFP';

%% Experiment info
depth = 9000; % um 
dz = 20; % spacing of channels (um)
gain = 50;
nChan = 384; % Number of channels
fs = 1500; % LFP sampling rate

%% Load data
lfp_file = dir(fullfile(lfp_dir,'*LFP_nt*.dat')); % list of LFP .dat files

lfp = NaN(960,fs*60); % sample only the initial 1 min
for ch=1:nChan
    str = extract(lfp_file(ch).name,digitsPattern);
    ch_depth = 10*(str2double(str{end-1})-1000);
    lfp_data = readTrodesExtractedDataFile(fullfile(lfp_file(ch).folder, lfp_file(ch).name)); % Load LFP data
    lfp(ch_depth/10,:) = double(lfp_data.fields.data(1:fs*60));
end

% Subtract channel and common noises
lfp = lfp - median(lfp,2,'omitnan');
lfp = lfp - median(lfp,1,'omitnan');
% lfp = lfp / gain; % convert to µV


figure
imagesc(corr(double(lfp')))
title('LFP Correlation')
xlabel('Channel ID')
ylabel('Channel ID')
% yticklabels(yticks*20)
colormap jet