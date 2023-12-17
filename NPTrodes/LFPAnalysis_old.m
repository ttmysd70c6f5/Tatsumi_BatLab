% function LFP_Analysis_v2()
%% Select directory
lfp_dir = 'C:\Users\Tatsumi\Documents\KQTY_NP\FlyingPixels\29968\230531\ephys\20230531_115208_merged.LFP';

%% Experiment info
depth = 9000; % um 
dz = 20; % spacing of channels (um)
gain = 50;

%% Load data
nChan = 384; % Number of channels
fs = 1500; % LFP sampling rate
fs_sp = 30000; % spikes sampling rate

lfp_file = dir(fullfile(lfp_dir,'*LFP_nt*.dat')); % list of LFP .dat files
tt = dir(fullfile(lfp_dir,'*timestamps.dat'));
tt = readTrodesExtractedDataFile(fullfile(tt.folder,tt.name)); % Load timestamps for LFP
tt = double(tt.fields.data); % timestamp (point)
tt = tt - tt(1) + fs_sp/fs; % offset timestamps
tt = tt / fs_sp;
rec_len = length(tt); % Recording time (samples)

startup

lfp = zeros(nChan,rec_len);
for chMaxChIdx=1:nChan
    lfp_data = readTrodesExtractedDataFile(fullfile(lfp_file(MaxChIdx).folder, lfp_file(MaxChIdx).name)); % Load LFP data
    lfp(MaxChIdx,:) = double(lfp_data.fields.data);
end

% Subtract channel and common noises
lfp = lfp - median(lfp,2);
lfp = lfp - median(lfp,1);
lfp = lfp / gain; % convert to µV

% % Sample initial 1 min
lfp = lfp(:,1:fs*60);
tt = tt(1:fs*60);

RecDate = strsplit(lfp_file(1).folder,{'\'});
RecDate = strsplit(RecDate{end},{'_'});
RecDate = RecDate{1};

%% LFP correlation
figure
set(gcf,'units','normalized','position',[0.1 0.1 0.7 0.5])
tl = tiledlayout(1,2,'TileSpacing','tight','Padding','compact');

% no filter
nexttile
imagesc(corr(double(lfp')))
title('LFP Correlation (raw)')
xlabel('Channel ID')
ylabel('From tip (mm)')
yticklabels(yticks*20)
colormap jet

% % remove highly deviated samples
% lfp_corrected = lfp - smoothdata(lfp,2);
% 
% figure; hold on
% plot(lfp(1:30:384,:)','k')
% plot(smoothdata(lfp(1:30:384,:),2)','y')
% plot(lfp_corrected(1:30:384,:)','b')
% plot(std(lfp_corrected(1:30:384,:)),'r')
% hold off
% % figure;plot(lfp_filtered(1:30:384,:)')
% 
% lfp_corrected = lfp_corrected(:,std(lfp_corrected(1:30:384,:))<4);
% 
% figure
% imagesc(corr(double(lfp_corrected')))
% title('LFP Correlation (high deviated points excluded)')
% xlabel('Channel ID')
% ylabel('From tip (mm)')
% yticklabels(yticks*20)
% colormap jet

% notch filter
x = lfp;
Fs = 1500;
d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',Fs);

% fvtool(d,'Fs',Fs)
x_filtered = filtfilt(d,x);

nexttile
imagesc(corr(double(x_filtered')))
title('LFP Correlation (notch filtered)')
xlabel('Channel ID')
ylabel('From tip (mm)')
yticklabels(yticks*20)
colormap jet

%% LFP correlation
%===FIGURE: probe correlation
figure;
set(gcf,'units','normalized','position',[0.01 0.55 0.25 0.35])

% nexttile(1)

% imagesc(corr(double(lfp')))
imagesc(corr(double(lfp')))
title('LFP Correlation')
xlabel('Channel ID')
ylabel('From tip (mm)')
yticklabels(yticks*20)
colormap jet

%% Current source density
% % estimate CSD
% % 5 point hamming filter from Ulbert et al. J Neurosci Methoods 2001 
% % ('Multiple microelectrode-recording system for intracortical 
% % applications') - equation 5 for spatial smoothing of signal
% clear CSD
% w = [0.23; 0.08; -0.62; 0.08; .23];
% for i = 1:size(lfp,1)
%     if i-2>0 && i+2<size(lfp,1)+1
%         u1 = lfp(i-2,:);
%         u2 = lfp(i-1,:);
%         u3 = lfp(i,:);
%         u4 = lfp(i+1,:);
%         u5 = lfp(i+2,:);
%         
%         CSD(i,:) = -(w(1)*u1 + w(2)*u2 + w(3)*u3 + w(4)*u4 +w(5)*u5)/(2*dz*2*dz);
%     end
% end
% CSD = CSD(3:end,:);
% 
% %plot LFP and CSD as function of depth
% % z = [depth-dz*(nChan-1):dz:depth];
% z = dz:dz:dz*(nChan-1);
% v=1; % interpolation factor
% 
% if isnan(z)
%     z = [1:1:nChan];
% end
% 
% % figure
% % subplot  121
% % imagesc(tt, z, interp2(lfp,v))
% % ylim([z(3) z(end-2)])
% % xlabel('time (s)')
% % ylabel('depth')
% % title('LFP (\muV)')
% 
% % subplot 122
% nexttile(2,[1 2])
% imagesc(tt, z(3:end-2), interp2(CSD,v))
% % imagesc(interp2(CSD,v))
% xlabel('time (s)')
% ylabel('From tip (µm)')
% colormap('jet')
% title('CSD')
% 
% h = colorbar;
% h.Ticks =  h.Limits;
% h.TickLabels = {'sink' 'source'};

%% Save
mkdir(fullfile(pwd,'Analysis'))
% saveas(tl,fullfile(pwd,'Analysis','lfp_summary_1'),'png')
exportgraphics(gcf,fullfile(pwd,'Analysis','lfp_summary_1.png'))


%% FIGURE2 setting
figure
set(gcf,'units','normalized','position',[0.01 0.1 0.95 0.8])
tl = tiledlayout(3,5,'TileSpacing','tight','Padding','compact');

%% LFP power band
% %===FIGURE: Delta (2-4Hz)
% lfp_delta = bandpass(lfp',[2 4],fs)'; % Extract high-frequency ripples
% delta_power = abs(hilbert(lfp_delta')'); % Hilbert transform of lfp

% Theta (4-8Hz)
lfp_theta = bandpass(lfp',[4 8],fs)'; % Extract high-frequency ripples
theta_power = abs(hilbert(lfp_theta')'); % Hilbert transform of lfp
theta_pwmean = mean(theta_power,2);

% Gamma (30-100Hz)
lfp_gamma = bandpass(lfp',[30 100],fs)'; % Extract high-frequency ripples
gamma_power = abs(hilbert(lfp_gamma')'); % Hilbert transform of lfp
gamma_pwmean = mean(gamma_power,2);

% Ripple (80-160Hz)
lfp_ripple = bandpass(lfp',[80 160],fs)'; % Extract high-frequency ripples
ripple_power = abs(hilbert(lfp_ripple')'); % Hilbert transform of lfp
ripple_pwmean = mean(ripple_power,2);

%===FIGURE
% figure
% set(gcf,'units','normalized','position',[0.265 0.1 0.45 0.8])
% tl = tiledlayout(3,2,'TileSpacing','tight','Padding','compact');

% Theta
nexttile(1)
imagesc(tt,1:nChan,log10(theta_power))
colormap jet
colorbar
% ylabel('Channel ID')
ylabel('From tip (µm)')
yticklabels(yticks*20)
xlabel('Time (sec)')
title('Log10 theta power (4-8Hz)')

nexttile(2)
plot(movmean(log10(theta_pwmean),5),'k')
title('Theta power band')
% xlabel('Channel ID')
xlabel('From tip (µm)')
xticklabels(xticks*20)
ylabel('Log10 mean LFP power (µV^2)')

% Gamma
nexttile(6)
imagesc(tt,1:nChan,log10(gamma_power))
colormap jet
colorbar
% ylabel('Channel ID')
ylabel('From tip (µm)')
yticklabels(yticks*20)
xlabel('Time (sec)')
title('Log10 gamma power (30-100Hz)')

nexttile(7)
plot(movmean(log10(gamma_pwmean),5),'k')
title('Gamma power band')
% xlabel('Channel ID')
xlabel('From tip (µm)')
xticklabels(xticks*20)
ylabel('Log10 mean LFP power (µV^2)')

% Ripple
nexttile(11)
imagesc(tt,1:nChan,log10(ripple_power))
colormap jet
colorbar
% ylabel('Channel ID')
ylabel('From tip (µm)')
yticklabels(yticks*20)
xlabel('Time (sec)')
title('Log10 ripple power (80-160Hz)')

nexttile(12)
plot(movmean(log10(ripple_pwmean),5),'k')
title('Ripple power band')
% xlabel('Channel ID')
xlabel('From tip (µm)')
xticklabels(xticks*20)
ylabel('Log10 mean LFP power (µV^2)')

%% Theta peaks
% Peaks in gamma band
th_thpks = 0.1;
[pks,locs,~,prominence] = findpeaks(movmean(log10(theta_pwmean),10),'MinPeakProminence',th_thpks);
% disp('Peaks and locs are')
% disp([pks,locs])

theta_mvmean = movmean(log10(theta_pwmean),5);
% pks_width = zeros(2,length(pks)); % index for width of peaks. 1st row is start and 2nd row is end of peaks
pks_width = struct([]);
for i = 1:length(pks)
    thetapks = theta_mvmean > pks(i)/2;
    [pks_start,pks_end] = findchanges(thetapks,0,1);
    pks_start(locs(i):end) = false; % Remove the start of peaks after the targe peak
    pks_end(1:locs(i)) = false; % Remove the end of peaks before the targe peak
    pks_start = find(pks_start); pks_end = find(pks_end);
    if isempty(pks_end)
        pks_end = nChan; % If the peak is at the right edge
    end
    if isempty(pks_start)
        pks_end = 1; % If the peak is at the left edge
    end

    [~,closestIdx] = min(abs(locs(i)-pks_start)); % Find start index of i-th peak
%     pks_width(1,i) = pks_start(closestIdx);
    pks_width(i).start = pks_start(closestIdx);
    [~,closestIdx] = min(abs(locs(i)-pks_end)); % Find start index of i-th peak
%     pks_width(2,i) = pks_end(closestIdx);
    pks_width(i).end = pks_end(closestIdx);
end

% figure
% set(gcf,'units','normalized','position',[0.74 0.63 0.21 0.25])
nexttile(3)

% figure
plot(movmean(log10(theta_pwmean),5),'k')
hold on
for i=1:length(pks_width)
    plot([pks_width(i).start,pks_width(i).end],[pks(i)/2 pks(i)/2],'--b')
    plot([locs(i) locs(i)],[pks(i) pks(i)-prominence(i)],'--b')
end
hold off
title('Theta peaks')
xlabel('From tip (µm)')
xticklabels(xticks*20)
ylabel('Log10 mean LFP power (µV^2)')

%% Theta power spectrum
params.tapers=[5 9];
params.pad=0;
params.Fs=1500;
params.fpass=[0 300];
params.trialave=0;

[S,f] = mtspectrumc(lfp',params); % Chronux toolbox

%===FIGURE
% figure
% set(gcf,'units','normalized','position',[0.74 0.05 0.21 0.5])
% tl = tiledlayout(length(pks_width),1,'TileSpacing','tight','Padding','compact');

if ~isempty(pks_width)
    for i=1:length(pks_width)
        if i<=2
            nexttile(3+i)
            plot(f,movmean(mean(S(:,pks_width(i).start:pks_width(i).end),2),100),'k'); xlim([0 200])
            title(sprintf('From tip %d - %d (µm)',dz*pks_width(i).start,dz*pks_width(i).end))
            xlabel('Frequency (Hz)')
            ylabel('Power')
        end
    end
end

%% Gamma peaks
% Peaks in gamma band
th_gmpks = 0.1;
[pks,locs,~,prominence] = findpeaks(movmean(log10(gamma_pwmean),10),'MinPeakProminence',th_gmpks);
% disp('Peaks and locs are')
% disp([pks,locs])

gamma_mvmean = movmean(log10(gamma_pwmean),5);
% pks_width = zeros(2,length(pks)); % index for width of peaks. 1st row is start and 2nd row is end of peaks
pks_width = struct([]);
for i = 1:length(pks)
    gammapks = gamma_mvmean > pks(i)/2;
    [pks_start,pks_end] = findchanges(gammapks,0,1);
    pks_start(locs(i):end) = false; % Remove the start of peaks after the targe peak
    pks_end(1:locs(i)) = false; % Remove the end of peaks before the targe peak
    pks_start = find(pks_start); pks_end = find(pks_end);
    if isempty(pks_end)
        pks_end = nChan; % If the peak is at the right edge
    end
    if isempty(pks_start)
        pks_end = 1; % If the peak is at the left edge
    end

    [~,closestIdx] = min(abs(locs(i)-pks_start)); % Find start index of i-th peak
%     pks_width(1,i) = pks_start(closestIdx);
    pks_width(i).start = pks_start(closestIdx);
    [~,closestIdx] = min(abs(locs(i)-pks_end)); % Find start index of i-th peak
%     pks_width(2,i) = pks_end(closestIdx);
    pks_width(i).end = pks_end(closestIdx);
end

% figure
% set(gcf,'units','normalized','position',[0.74 0.63 0.21 0.25])
nexttile(8)

% figure
plot(movmean(log10(gamma_pwmean),5),'k')
hold on
for i=1:length(pks_width)
    plot([pks_width(i).start,pks_width(i).end],[pks(i)/2 pks(i)/2],'--b')
    plot([locs(i) locs(i)],[pks(i) pks(i)-prominence(i)],'--b')
end
hold off
title('Gamma peaks')
xlabel('From tip (µm)')
xticklabels(xticks*20)
ylabel('Log10 mean LFP power (µV^2)')

%% Gamma power spectrum
params.tapers=[5 9];
params.pad=0;
params.Fs=1500;
params.fpass=[0 300];
params.trialave=0;

[S,f] = mtspectrumc(lfp',params); % Chronux toolbox

%===FIGURE
% figure
% set(gcf,'units','normalized','position',[0.74 0.05 0.21 0.5])
% tl = tiledlayout(length(pks_width),1,'TileSpacing','tight','Padding','compact');

if ~isempty(pks_width)
    for i=1:length(pks_width)
        if i<=2
            nexttile(8+i)
            plot(f,movmean(mean(S(:,pks_width(i).start:pks_width(i).end),2),100),'k'); xlim([0 200])
            title(sprintf('From tip %d - %d (µm)',dz*pks_width(i).start,dz*pks_width(i).end))
            xlabel('Frequency (Hz)')
            ylabel('Power')
        end
    end
end

%% Ripple peaks
% Peaks in gamma band
th_rplpks = 0.1;
[pks,locs,~,prominence] = findpeaks(movmean(log10(ripple_pwmean),10),'MinPeakProminence',th_rplpks);
% disp('Peaks and locs are')
% disp([pks,locs])

ripple_mvmean = movmean(log10(ripple_pwmean),5);
% pks_width = zeros(2,length(pks)); % index for width of peaks. 1st row is start and 2nd row is end of peaks
pks_width = struct([]);
for i = 1:length(pks)
    ripplepks = ripple_mvmean > pks(i)/2;
    [pks_start,pks_end] = findchanges(ripplepks,0,1);
    pks_start(locs(i):end) = false; % Remove the start of peaks after the targe peak
    pks_end(1:locs(i)) = false; % Remove the end of peaks before the targe peak
    pks_start = find(pks_start); pks_end = find(pks_end);
    if ~isempty(pks_start) && isempty(pks_end)
        pks_end = nChan; % If the peak is at the right edge
    elseif isempty(pks_start) && ~isempty(pks_end)
        pks_start = 1; % If the peak is at the left edge
    end

    if ~isempty(pks_start) && ~isempty(pks_end)
        [~,closestIdx] = min(abs(locs(i)-pks_start)); % Find start index of i-th peak
    %     pks_width(1,i) = pks_start(closestIdx);
        pks_width(i).start = pks_start(closestIdx);
        [~,closestIdx] = min(abs(locs(i)-pks_end)); % Find start index of i-th peak
    %     pks_width(2,i) = pks_end(closestIdx);
        pks_width(i).end = pks_end(closestIdx);
    end
end

% figure
% set(gcf,'units','normalized','position',[0.74 0.63 0.21 0.25])
nexttile(13)

% figure
plot(movmean(log10(ripple_pwmean),5),'k')
hold on
for i=1:length(pks_width)
    plot([pks_width(i).start,pks_width(i).end],[pks(i)/2 pks(i)/2],'--b')
    plot([locs(i) locs(i)],[pks(i) pks(i)-prominence(i)],'--b')
end
hold off
title('Ripple peaks')
xlabel('From tip (µm)')
xticklabels(xticks*20)
ylabel('Log10 mean LFP power (µV^2)')

%% Ripple power spectrum
params.tapers=[5 9];
params.pad=0;
params.Fs=1500;
params.fpass=[0 300];
params.trialave=0;

[S,f] = mtspectrumc(lfp',params); % Chronux toolbox

%===FIGURE
% figure
% set(gcf,'units','normalized','position',[0.74 0.05 0.21 0.5])
% tl = tiledlayout(length(pks_width),1,'TileSpacing','tight','Padding','compact');

if ~isempty(pks_width)
    for i=1:length(pks_width)
        if i<=2
            nexttile(13+i)
            plot(f,movmean(mean(S(:,pks_width(i).start:pks_width(i).end),2),100),'k'); xlim([0 200])
            title(sprintf('From tip %d - %d (µm)',dz*pks_width(i).start,dz*pks_width(i).end))
            xlabel('Frequency (Hz)')
            ylabel('Power')
        end
    end
end


%% Save
mkdir(fullfile(pwd,'Analysis'))
% saveas(tl,fullfile(pwd,'Analysis','lfp_summary_1'),'png')
exportgraphics(gcf,fullfile(pwd,'Analysis','lfp_summary_2.png'))


%% Putative ripples
% Ripple (80-160Hz) detection
lfp_ripple = bandpass(lfp',[80 160],fs)'; % Extract high-frequency ripples
ripple_power = abs(hilbert(lfp_ripple')'); % Hilbert transform of lfp
ripple_pwmean = mean(ripple_power,2); % Mean of lfp power
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




%% Save
mkdir(fullfile(pwd,'Analysis'))
% saveas(tl,fullfile(pwd,'Analysis','lfp_summary_1'),'png')
exportgraphics(gcf,fullfile(pwd,'Analysis','lfp_summary_3.png'))

%% Granger causality
% y1 = lfp_gamma(150,:)';
% y2 = lfp_gamma(350,:)';
% 
% gctest_kpss = kpsstest(y1) & kpsstest(y2);
% [h,p,stat,cvalue] = gctest(y1,y2,NumLags=1,Alpha=0.01);
% [h2,p2,stat2,cvalue2] = gctest(y2,y1,NumLags=1,Alpha=0.01);
% 
% % save(fullfile(cd,'gctest.mat'),'gctest_h','gctest_p','gctest_stat','gctest_c','gctest_kpss')

%% Skip noisy channels
% lfpcorr = corr(double(lfp'));
% % corrmean = mean(lfpcorr,2);
% corrstd = std(lfpcorr,0,2);
% % corrstd = mad(lfpcorr,1,2);
% corrdiff = diff(corrstd);
% % corrdiffsd = std(corrdiff);
% corrdiffsd = 2*mad(corrdiff,1);
% % corrdiffmean = mean(corrdiff);
% corrdiffmean = median(corrdiff);
% 
% 
% figure
% set(gcf,'units','normalized','position',[0.35 0.1 0.5 0.3])
% tiledlayout(1,2,'TileSpacing','tight','Padding','compact');
% 
% %===FIGURE: detection of noisy channels
% nexttile
% plot(corrdiff)
% hold on
% plot(xlim,[corrdiffmean+corrdiffsd,corrdiffmean+corrdiffsd],'r')
% plot(xlim,[corrdiffmean-corrdiffsd,corrdiffmean-corrdiffsd],'r')
% hold off
% 
% %===FIGURE: std of LFP signals
% steps = abs(corrdiff - corrdiffmean) > corrdiffsd;
% corrdiff_ip = interp1(find(~steps),corrstd(~steps),1:length(corrstd));
% 
% nexttile
% plot(corrstd)
% hold on
% plot(corrdiff_ip)
% hold off

%% LFP power with different band
% 
% % Sample the first 1 min
% lfp = lfp(:,1:fs*60);
% tt = tt(1:fs*60);
% 
% th = 7;
% 
% figure
% set(gcf,'units','normalized','position',[0.01 0.1 0.98 0.3])
% tl = tiledlayout(1,5,'TileSpacing','tight','Padding','compact');
% 
% %===FIGURE: probe correlation
% ax(1) = nexttile;
% imagesc(corr(double(lfp')))
% title('Probe Correlation')
% xlabel('Depth (mm)')
% ylabel('Channel ID')
% xticklabels(xticks*2)
% colormap jet
% 
% 
% %===FIGURE: Theta (4-8Hz)
% lfp_theta = bandpass(lfp',[4 8],fs)'; % Extract high-frequency ripples
% lfp_power = abs(hilbert(lfp_theta')'); % Hilbert transform of lfp
% lfp_pwmean = mean(lfp_power,2); % Mean of lfp power
% lfp_pwsd = std(lfp_power,0,2);
% 
% Theta = false(size(lfp_power));
% for ch=1:nChan
%     Theta(ch,:) = lfp_power(ch,:) > lfp_pwmean(ch) + th*lfp_pwsd(ch);
% end
% 
% nexttile
% [theta_init,theta_end] = findchanges(double(Theta),0,2);
% plot(sum(theta_init,2),'k');
% title('Theta: power (4-8Hz) > 7SD')
% xlabel('Channel ID')
% ylabel('# of events')
% xticklabels(xticks*2)
% xlabel('Depth (mm)')
% % ylabel('# of putative ripples')
% xlim([1,nChan])
% 
% %===FIGURE: Delta (2-4Hz)
% lfp_delta = bandpass(lfp',[2 4],fs)'; % Extract high-frequency ripples
% lfp_power = abs(hilbert(lfp_delta')'); % Hilbert transform of lfp
% lfp_pwmean = mean(lfp_power,2); % Mean of lfp power
% lfp_pwsd = std(lfp_power,0,2);
% 
% Delta = false(size(lfp_power));
% for ch=1:nChan
%     Delta(ch,:) = lfp_power(ch,:) > lfp_pwmean(ch) + 7*lfp_pwsd(ch);
% end
% 
% nexttile
% [delta_init,delta_end] = findchanges(double(Delta),0,2);
% plot(sum(delta_init,2),'k');
% title('Delta: power (2-4Hz) > 7SD')
% xlabel('Channel ID')
% xticklabels(xticks*2)
% xlabel('Depth (mm)')
% xlim([1,nChan])
% 
% 
% %===FIGURE: Gamma (30-100Hz)
% lfp_gamma = bandpass(lfp',[30 100],fs)'; % Extract high-frequency ripples
% lfp_power = abs(hilbert(lfp_gamma')'); % Hilbert transform of lfp
% lfp_pwmean = mean(lfp_power,2); % Mean of lfp power
% lfp_pwsd = std(lfp_power,0,2);
% 
% Gamma = false(size(lfp_power));
% for ch=1:nChan
%     Gamma(ch,:) = lfp_power(ch,:) > lfp_pwmean(ch) + th*lfp_pwsd(ch);
% end
% 
% nexttile
% [gamma_init,gamma_end] = findchanges(double(Gamma),0,2);
% gamma_cnt = sum(gamma_init,2);
% % steps = abs(corrdiff - corrdiffmean) > corrdiffsd;
% % gamma_cnt = interp1(find(~steps),gamma_cnt(~steps),1:length(corrstd)); % remove abrupt changes
% 
% plot(gamma_cnt,'k');
% title('Gamma: power (30-100Hz) > 7SD')
% xlabel('Channel ID')
% xticklabels(xticks*2)
% xlabel('Depth (mm)')
% xlim([1,nChan])
% 
% 
% %===FIGURE: # of ripples
% lfp_filtered = bandpass(lfp',[80 160],fs)'; % Extract high-frequency ripples
% lfp_power = abs(hilbert(lfp_filtered')'); % Hilbert transform of lfp
% lfp_pwmean = mean(lfp_power,2); % Mean of lfp power
% lfp_pwsd = std(lfp_power,0,2);
% 
% Ripples = false(size(lfp_power));
% for ch=1:nChan
%     Ripples(ch,:) = lfp_power(ch,:) > lfp_pwmean(ch) + th*lfp_pwsd(ch);
% end
% 
% nexttile
% [ripples_init,ripples_end] = findchanges(double(Ripples),0,2);
% ripples_cnt = sum(ripples_init,2);
% % steps = abs(corrdiff - corrdiffmean) > corrdiffsd;
% % ripples_cnt = interp1(find(~steps),ripples_cnt(~steps),1:length(corrstd)); % remove abrupt changes
% 
% plot(ripples_cnt,'k');
% title('Ripples: power (80-160Hz) > 7SD')
% xlabel('Channel ID')
% xticklabels(xticks*2)
% xlabel('Depth (mm)')
% xlim([1,nChan])



%% LFP coherence between putative MEC and CA1

%% LFP traces

%% Putative ripples
% lfp_power = abs(hilbert(lfp_filtered')'); % Hilbert transform of lfp
% lfp_pwmean = mean(lfp_power,2); % Mean of lfp power
% lfp_pwsd = std(lfp_power,0,2);
% 
% Ripples = false(size(lfp_power));
% for ch=1:nChan
%     Ripples(ch,:) = lfp_power(ch,:) > lfp_pwmean(ch) + 7*lfp_pwsd(ch);
% end
% 
%  
% figure
% set(gcf,'units','normalized','position',[0.35 0.05 0.30 0.8])
% tl = tiledlayout(2,1,'TileSpacing','tight','Padding','compact');
% %===FIGURE: probe correlation
% ax(1) = nexttile;
% imagesc(corr(double(lfp')))
% title('Probe Correlation')
% xlabel('Channel ID')
% ylabel('Channel ID')
%  
%  
% %===FIGURE: # of ripples
% [ripples_init,ripples_end] = findchanges(double(Ripples),0,2);
% ax(2) = nexttile;
% plot(sum(ripples_init,2),'k');
% title('LFP power (80-160Hz) > 7SD')
% xlabel('Channel ID')
% ylabel('# of putative ripples')
% xlim([1,nChan])
% % linkaxes(ax)
% 
% %===FIGURE: Example of putative ripples
% ch=240;
% idx = find(ripples_init(ch,:));
% Nplots = 5;
% idx_extracted = idx(randperm(length(idx),Nplots));
% figure
% set(gcf,'units','normalized','position',[0.05 0.05 0.9 0.8])
% tl = tiledlayout(Nplots,2);
% title(tl,sprintf('Channel %d',ch))
% xlabel(tl,'Time (s)')
% ylabel(tl,'Amplitude')
% for i=1:Nplots
%     if idx_extracted(i) > fs*10 && idx_extracted(i) < size(lfp_filtered,2) - fs*10
%         nexttile
%         plot(lfp_filtered(ch,idx_extracted(i)-fs*10:idx_extracted(i)+fs*10),'k')
%         xticks(1:fs*10:fs*20+1)
%         xticklabels((xticks-1)/fs-10)
%         xlim([1,fs*20+1])
%         title('Filtered LFP (80-160Hz)')
%         
%         nexttile
%         plot(lfp(ch,idx_extracted(i)-fs*10:idx_extracted(i)+fs*10),'k')
%         xticks(1:fs*10:fs*20+1)
%         xticklabels((xticks-1)/fs-10)
%         xlim([1,fs*20+1])
%         title('Unfiltered LFP')
%     end
% end
%  
% %===FIGURE: average waveform
% wv = zeros(length(idx),fs*0.2+1);
% for i = 1:length(idx)
%     wv(i,:) = lfp_filtered(ch,idx(i)-fs*0.1:idx(i)+fs*0.1);
% end
% wv_mean = mean(wv,1);
% 
% figure
% plot(-100:1000*1/fs:100,wv_mean,'k')
% xlabel('Time (ms)')
% ylabel('Amplitude')




%% Spectrogram of LFP
% Nfft = 2^nextpow2(length(lfp_filtered));
% freqfft = (0:(Nfft/2-1)*(fs/Nfft));
% lfp_fft = fft(lfp_filtered,Nfft);
% 
% Nspec = 256;
% wspec = hamming(Nspec);
% Noverlap = Nspec/2;
% % Noverlap = Nspec-1;
% 
% [spect_lfp,fspec,tspec] = spectrogram(lfp_fft,wspec,Noverlap,Nspec,fs);
% spectrogram(lfp_fft,wspec,Noverlap,Nspec,fs,'yaxis');