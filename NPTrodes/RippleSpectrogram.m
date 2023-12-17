%% Experiment info
gain = 50; % LFP gain
nChan = 384; % Number of channels
fs = 1500; % LFP sampling rate
fs_sp = 30000; % spikes sampling rate



%% Select lfp directory
% lfp_dir = 'C:\Users\YartsevPixels1\Desktop\Tatsumi\Data\230531\20230531_115208_merged.LFP';
lfp_dir = 'C:\Users\Tatsumi\Documents\KQTY_NP\FlyingPixels\29968\230531\ephys\20230531_115208_merged.LFP';


%% Load data
lfp_file = dir(fullfile(lfp_dir,'*LFP_nt*.dat')); % list of LFP .dat files
tt = dir(fullfile(lfp_dir,'*timestamps.dat'));
tt = readTrodesExtractedDataFile(fullfile(tt.folder,tt.name)); % Load timestamps for LFP
tt = double(tt.fields.data); % timestamp (point)
tt = tt - tt(1) + fs_sp/fs; % offset timestamps
tt = tt / fs_sp;
rec_len = length(tt); % Recording time (samples)

pat = digitsPattern;
str = lfp_file(1).name;
ch_depth = zeros(nChan,1);

lfp = zeros(nChan,rec_len);
for ch=1:nChan
    lfp_data = readTrodesExtractedDataFile(fullfile(lfp_file(ch).folder, lfp_file(ch).name)); % Load LFP data
    lfp(ch,:) = double(lfp_data.fields.data);
    str = extract(lfp_file(ch).name,digitsPattern);
    ch_depth(ch,1) = 10*(str2double(str{3})-1000);
end


% Subtract channel and common noises
lfp = lfp - median(lfp,2);
lfp = lfp - median(lfp,1);
lfp = lfp / gain; % convert to µV

RecDate = strsplit(lfp_file(1).folder,{'\'});
RecDate = strsplit(RecDate{end},{'_'});
RecDate = RecDate{1};


%% FIGURE
%% temporal LFP power traces
% figure
% imagesc(tt,1:nChan,log10(ripple_power))
% colormap jet
% colorbar
% ylabel('Channel ID')
% xlabel('Time (sec)')
% title('Log10 ripple power (80-160Hz)')
% 
% %% Mean power
% figure
% plot(movmean(log10(ripple_pwmean),5),'k')
% title('Ripple power band')
% xlabel('Channel ID')
% ylabel('Log10 mean LFP power (µV^2)')

%% Spectrogram
ch = 1;

params.tapers=[5 9];
params.pad=0;
params.Fs=1500;
params.fpass=[0 300];
params.trialave=0;

svdir = 'C:\Users\YartsevPixels1\Desktop\Tatsumi\Data\230531\20230531_115208_merged.LFP\ripples';

nChan_page = 24;
nPages = ceil(nChan/nChan_page);

for page = 1:nPages
    f = figure;
    f.PaperOrientation = 'landscape';
    f.PaperPosition = [0.05 0.05 0.9 0.9];
    f.PaperUnits = 'normalized';
    set(gcf,'Units','normalized','OuterPosition',[0.05 0.05 0.9 0.9])
    tl = tiledlayout(4,6,'TileSpacing','Compact');

    ch_now = 0;
    while ch_now < nChan_page
        [S,f] = mtspectrumc(lfp(ch,:),params); % Chronux toolbox
    
        nexttile
        plot(f,log10(movmean(S,1000)),'k'); 
        xlim([0 200])
%         ylim([-0.1 0.5])
        title(sprintf('Electrode %d',ch_depth(ch)/10))
        xlabel('Frequency (Hz)')
        ylabel('Log10 power')

        pause(0.3)

        ch = ch + 1;
        ch_now = ch_now + 1;
    end

%     saveas(f,[svdir,'\ripples_',num2str(page+1000),'.pdf'])
    exportgraphics(gcf,[svdir,'\ripples_',num2str(page+1000),'.pdf'])
    pause(0.3)
    close all
end
%% Putative ripples
% Extract ripple (80-160Hz)
lfp_ripple = bandpass(lfp',[80 160],fs)'; % Extract high-frequency ripples
ripple_power = abs(hilbert(lfp_ripple')'); % Hilbert transform of lfp
ripple_pwmean = mean(ripple_power,2);
ripple_pwsd = std(ripple_power,0,2); % std of lfp power

th = 7;

Ripples = false(size(ripple_power));
for ch = 1:nChan
    Ripples(ch,:) = ripple_power(ch,:) > ripple_pwmean(ch) + th*ripple_pwsd(ch);
end

[~,MaxChIdx] = max(movmean(log10(ripple_pwmean),5));
[Ripples_ch,~] = findchanges(Ripples(MaxChIdx,:),0,2);
Ripples_idx = find(Ripples_ch);

ripple_wv = zeros(1,fs*0.2+1);

%%
% figure
% set(gcf,'units','normalized','position',[0.74 0.05 0.21 0.5])
% tl = tiledlayout(10,ceil(length(Ripples_idx)/10),'TileSpacing','tight','Padding','compact');

idx_pre = 0; % initialize parameter
ripple_cnt = 0; % # of ripples
for i=1:length(Ripples_idx)
    idx_now = Ripples_idx(i);
    
    if idx_now - idx_pre > fs*0.1 && idx_now + idx_pre < size(Ripples,2) - fs*0.1
%         nexttile
%         plot(-fs*0.1:fs*0.1, lfp_ripple(ch,Ripples_idx(i)-fs*0.1:Ripples_idx(i)+fs*0.1),'k')
%         title([num2str(Ripples_idx(i)/1500),' (sec)'])
    
        ripple_wv = ripple_wv + lfp_ripple(MaxChIdx,Ripples_idx(i)-fs*0.1:Ripples_idx(i)+fs*0.1);
        ripple_cnt = ripple_cnt + 1;
    end
    
    idx_pre = idx_now;
end

% figure
figure
plot(-fs*0.1:fs*0.1,ripple_wv,'k')
title(sprintf('putative ripple: %d µm from tip, %d counts',MaxChIdx*20,ripple_cnt))
ylabel('Amplitude')
xlabel('Time (ms)')
xlim([-0.1*fs 0.1*fs])


