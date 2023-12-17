function RippleSpectrum_LowMemory(ch,lfp_dir)
%% Experiment info
gain = 50; % LFP gain
fs = 1500; % LFP sampling rate

%% Load data
lfp_file = dir(fullfile(lfp_dir,'*LFP_nt*.dat')); % list of LFP .dat files
lfp_data = readTrodesExtractedDataFile(fullfile(lfp_file(ch).folder, lfp_file(ch).name)); % Load LFP data
lfp = double(lfp_data.fields.data);
str = extract(lfp_file(ch).name,digitsPattern);
ch_depth = 10*(str2double(str{3})-1000);

% Subtract median
load(fullfile(lfp_file(1).folder,'lfpmedian.mat'),'lfpmed1')
lfp = lfp - lfpmed1;
lfp = lfp - median(lfp);

% low pass
lfp = lowpass(lfp,400,fs);

% lfp = lfp / gain; % convert to µV

%%
params.tapers=[5 9];
params.pad=0;
params.Fs=fs;
params.fpass=[0 300];
params.trialave=0;

[S,f] = mtspectrumc(lfp,params); % Chronux toolbox

% figure
plot(f,log10(movmean(S,10000)),'k'); 
xlim([0 200])
% ylim([-0.1 0.5])
title(sprintf('Electrode %d',ch_depth/10))
xlabel('Frequency (Hz)')
ylabel('Log10 power')