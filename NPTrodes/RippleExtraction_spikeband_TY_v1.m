%% Function for the analysis of Extracellular Field Recordings in the Hippocampus and the detection of Ripples
% Original script is written by Angelo Forli (TT_CSC_Analyze_AF_v0.m)

lfpDir = 'C:\Users\Tatsumi\Documents\KQTY_NP\FlyingPixels\29968\230531\ephys\20230531_115208_merged.LFP';
spDir = 'C:\Users\Tatsumi\Documents\KQTY_NP\FlyingPixels\29968\230531\ephys\20230531_115208_merged.spikeband';
tsDir = 'C:\Users\Tatsumi\Documents\KQTY_NP\FlyingPixels\29968\230531\ephys\20230531_115208_merged.time';

lfpFile = dir(fullfile(lfpDir,'*LFP_nt*.dat')); % list of LFP .dat files
spFile = dir(fullfile(spDir,'*spikeband_nt*.dat')); % list of LFP .dat files

nChan = 384;

confirmed_by_inspection = zeros(nChan,1); % (# of electrodes)

%%
ch = 100;

fprintf('Analyzing electrode %d\n', ch)

%% Load lfp data
Fs = 1500;

lfpData = readTrodesExtractedDataFile(fullfile(lfpFile(ch).folder, lfpFile(ch).name)); % Load LFP data
lfp = double(lfpData.fields.data);
str = extract(lfpFile(ch).name,digitsPattern);
ch_depth = 10*(str2double(str{3})-1000);

% time train
tsFile = dir(fullfile(lfpDir,'*timestamps.dat'));
tstamps = readTrodesExtractedDataFile(fullfile(tsFile.folder,tsFile.name)); % Load timestamps for LFP
tstamps = double(tstamps.fields.data); % timestamp (point)
tstamps = tstamps - tstamps(1) + Fs_sp/Fs; % offset timestamps
tstamps = tstamps / Fs_sp;

% Subtract median
load(fullfile(lfpFile(1).folder,'lfpmedian.mat'),'lfpmed1')
lfp = lfp - lfpMed1;
lfp = lfp - median(lfp);

%% Load spikeband data
Fs = 30000;
ext_len = Fs*5;
ts_start = Fs*10;

spData = readTrodesExtractedDataFile(fullfile(spFile(ch).folder, spFile(ch).name)); % Load spikeband data
sp = double(spData.fields.data);
str = extract(spFile(ch).name,digitsPattern);
ch_depth = 10*(str2double(str{3})-1000);

tsFile = dir(fullfile(tsDir,'*timestamps.dat'));
tstamps = readTrodesExtractedDataFile(fullfile(tsFile.folder,tsFile.name)); % Load timestamps for spikeband
tstamps = tstamps.fields(1).data;

% Subtract median
load(fullfile(spFile(1).folder,'spmedian.mat'),'spMed1')
% sp = sp - spMed1;
% sp = sp - median(sp);

sp = sp(ts_start+1:ts_start+ext_len) - spMed1;
sp = sp - median(sp);

%% Options
CSC_options.lookATpre = 0;
CSC_options.lookATfly = 0;
CSC_options.lookATpst = 0;
CSC_options.lookATsorted = 0;
CSC_options.checkAlgn = 1;
CSC_options.extractRipples = 1;
CSC_options.loadBehavior = 1;
CSC_options.lookAccel = 0;
CSC_options.savedata = 1;
CSC_options.savefigures = 0;  fig_count = 1;

RP_F_band = [100 200];                  % Ripples frequency band (From MY: [80 160])
RP_th = 5;                              % Threshold on zscored Ripple power (3.5)
dec = 10;                               % Downsampling factor for Ripple detection

%% Ripple detection
    
%=== Downsample raw voltage and timestamps
ds_Fs = Fs/dec;                                                     % Downsampled Frequency
ds_V = downsample(sp,dec);                                       % Downsampled Voltage
ds_tstamps = downsample(tstamps,dec);                               % Downsampled TimeStamps (s)

%=== Define filter for MU
SP_F_band = [600 6000];                 % Spikes frequency band
passband_norm = SP_F_band./(0.5*Fs);    % Normalize over sampling frequency
[b,a] = ellip(4,0.1,70,passband_norm);  % Define filter transfer function coefficients

    
%=== Ripple Detection (amplitude,time,duration)
ripple_signal = bandpass(ds_V,RP_F_band,Fs/dec);                    % Bandpass signal at the Ripple Band
ripple_power  = zscore(abs(hilbert(ripple_signal)));                % Calculate z-scored power as the magnitude of the Hilbert transform
ripple_power  = smoothdata(ripple_power,'gaussian',0.05*Fs/dec);    % Smooth with a 50 ms Gaussian kernel
[RPL_a,RPL_t,RPL_d] = findpeaks(ripple_power,Fs/dec,...             % Detect peaks in the ripple power...
                                'MinPeakHeight',RP_th,...           % larger than RP_th...
                                'MinPeakDistance',0.1);             % with minimal separation of 100 ms

%=== Remove all the ripples happening during flight (+-0.5s)
% extended_flight = movmax(bflying(:,imp_bat),[50 50]);
% r2k = true(size(RPL_t));
% parfor i=1:numel(RPL_t)
%     if RPL_t(i)>(TTL_timestamps.fly(1)-Timestamps_of_first_samples_usec(1))/1e6 && RPL_t(i)<(TTL_timestamps.fly(end)-Timestamps_of_first_samples_usec(1))/1e6
%         [~,idx] = min(abs((RPL_t(i)-(TTL_timestamps.fly(1)-Timestamps_of_first_samples_usec(1))/1e6)-t));
%         if extended_flight(idx,1)
%             r2k(i) = 0
%         end
%     end
% end
% RPL_a = RPL_a(r2k); RPL_t = RPL_t(r2k); RPL_d = RPL_d(r2k);
    
%=== Calculate range of the raw voltage (+-50ms around ripple peak) or maximum derivative
RPL_r = zeros(size(RPL_t));
for i=1:numel(RPL_t) 
    sample_ripple = round(RPL_t(i)*Fs);
    if sample_ripple-0.05*Fs>0 && sample_ripple+0.05*Fs<numel(raw_V)
        %RPL_r(1,i) = range(raw_V(round(sample_ripple-0.05*Fs):round(sample_ripple+0.05*Fs)));
        RPL_r(1,i) = max(abs(diff(raw_V(round(sample_ripple-0.05*Fs):round(sample_ripple+0.05*Fs)))))*Fs;
        RPL_r(1,i) = RPL_r(1,i)/1e6;    %uV/us
    end
end
    
%=== Remove some artifacts
RPL_t=RPL_t(RPL_r<5);   RPL_a=RPL_a(RPL_r<5); RPL_d=RPL_d(RPL_r<5); RPL_r=RPL_r(RPL_r<5);

%=== Plot a few examples (+-100ms around ripple peak)
[rs_RPL_t,rs_idx] = datasample(RPL_t,50);              % Select random subset
rs_RPL_a = RPL_a(rs_idx);
rs_RPL_d = RPL_d(rs_idx);
rs_RPL_r = RPL_r(rs_idx);
[~,st_idx] = sort(rs_RPL_r);                            % Sort according to...
rs_RPL_t = rs_RPL_t(st_idx);
rs_RPL_a = rs_RPL_a(st_idx);
rs_RPL_d = rs_RPL_d(st_idx);
rs_RPL_r = rs_RPL_r(st_idx);
figure('units','normalized','outerposition',[0 0 1 1]);
tiledlayout(5,20,'TileSpacing','tight');
for i =1:50
    sample_ripple = round(rs_RPL_t(i)*Fs);
    if round(sample_ripple-0.1*Fs)<0 || round(sample_ripple+0.1*Fs)>numel(raw_V),continue;end
    nexttile; 
    plot((round(-0.1*Fs):round(+0.1*Fs))/Fs,normalize(raw_V(round(sample_ripple-0.1*Fs):round(sample_ripple+0.1*Fs)),'range'),'k');
    hold on;    plot((round(-0.1*Fs):round(+0.1*Fs))/Fs,-1+normalize(filtfilt(b,a,raw_V(round(sample_ripple-0.1*Fs):round(sample_ripple+0.1*Fs))),'range')); hold off;
    xlim([-0.1 0.1]);   xlabel('s');    ylabel('uV');   %xticks([]); yticks([]);
    title(['A = ', num2str(rs_RPL_a(i),2), ' D = ' num2str(rs_RPL_d(i)*1e3,2),' ms, A = ', num2str(rs_RPL_r(i),2), ' uV/us']);
    nexttile;
    spectrogram_AF_v0(raw_V(round(sample_ripple-0.5*Fs):round(sample_ripple+0.5*Fs)),Fs,0.2,0.05,[70 200]);
    %plot(abs(diff(raw_V(round(sample_ripple-0.1*Fs):round(sample_ripple+0.1*Fs)))),'k');
    %     [cfs,f] = cwt(raw_V(round(sample_ripple-0.5*Fs):round(sample_ripple+0.5*Fs)),Fs);
    %     imagesc((round(-0.1*Fs):round(+0.1*Fs))/Fs,f(52:80),imgaussfilt(abs(cfs(52:80,:)),[2 500]));
    %     xlabel('Time (s)');ylabel('Frequency (Hz)');axis xy;
end
sgtitle(['TT ',num2str(TT),', CSC ',num2str(CSC_ch)]);
% fig_count = saveFig(analysis_directory,[batdate, '_CSC',  num2str(CSC_ch)],fig_count,CSC_options.savefigures);
    
%     %=== Plot whole session trace
%     figure('units','normalized','outerposition',[0 0.2 1 0.4]);
%     tiledlayout(2,1,'TileSpacing','tight');
%     ax(1) = nexttile;   plot((ds_tstamps-ds_tstamps(1))/60e6,ds_V); ylabel('Voltage (uV)');
%     text((TTL_timestamps.fly(1)-ds_tstamps(1))/60e6,2e3,'Start Fly Session');
%     text((TTL_timestamps.fly(end)-ds_tstamps(1))/60e6,2e3,'Stop Fly Session');
%     ax(2) = nexttile;   plot((ds_tstamps-ds_tstamps(1))/60e6,ripple_power); ylabel('Ripple Power (sd)');
%     xlabel('Minutes');  linkaxes(ax,'x');   sgtitle(['TT ',num2str(TT),', CSC ',num2str(CSC_ch)]);
%     fig_count = saveFig(analysis_directory,[batdate, '_CSC',  num2str(CSC_ch)],fig_count,CSC_options.savefigures);
    
    %% Plot spectrum associated with detected Ripples
    
    n=512; P=zeros(numel(RPL_t),n);
    for i = 1:numel(RPL_t)
        sample_ripple = round(RPL_t(i)*Fs/dec);
        if round(sample_ripple-0.5*Fs/dec)<0 || round(sample_ripple+0.5*Fs/dec)>numel(ds_V) || RPL_r(i)>5,continue;end  % Do not include artifacts and ripples at the edges of the session
        ripple_chunk = ds_V(sample_ripple-n/2:sample_ripple+n/2-1);
        Y = fft(ripple_chunk);
        P(i,:) = abs(Y/n).^2;
    end
    P_ave = mean(P,1);  f_ave = Fs/dec*(0:n/2)/n;
    % Y_all = fft(ds_V);
    % n_all = length(ds_V);
    % P_all = abs(Y_all/n_all).^2;
    % f_all = Fs/dec*(0:n_all/2)/n_all;
    
    figure;
    plot(f_ave,P_ave(1:n/2+1));
    xlim([50 200]);    set(gca, 'YScale', 'log');
    xlabel('Frequency (Hz)');   h = gca;    h.YAxis.Visible = 'off';
    sgtitle(['TT ',num2str(TT),', CSC ',num2str(CSC_ch)]);
    fig_count = saveFig(analysis_directory,[batdate, '_CSC',  num2str(CSC_ch)],fig_count,CSC_options.savefigures);
    
    %% Ask for Visual Inspection
    
    answer = questdlg('Are SWRs present on this Tetrode?');
    if strcmp(answer,'Yes')
        confirmed_by_inspection(TT) = 1;
    end
    
    %% Save detected Ripples, referenced to first TTL of the corresponding session
    
    idx_pre = find(RPL_t>(TTL_timestamps.pre(1)-Timestamps_of_first_samples_usec(1))/1e6 & RPL_t<(TTL_timestamps.pre(end)-Timestamps_of_first_samples_usec(1))/1e6);
    idx_fly = find(RPL_t>(TTL_timestamps.fly(1)-Timestamps_of_first_samples_usec(1))/1e6 & RPL_t<(TTL_timestamps.fly(end)-Timestamps_of_first_samples_usec(1))/1e6);
    idx_pst = find(RPL_t>(TTL_timestamps.pst(1)-Timestamps_of_first_samples_usec(1))/1e6 & RPL_t<(TTL_timestamps.pst(end)-Timestamps_of_first_samples_usec(1))/1e6);
    RPL_pre.t = RPL_t(idx_pre)'-(TTL_timestamps.pre(1)-Timestamps_of_first_samples_usec(1))/1e6;    RPL_pre.a = RPL_a(idx_pre)';    RPL_pre.d = RPL_d(idx_pre)';    RPL_pre.r = RPL_r(idx_pre)';
    RPL_fly.t = RPL_t(idx_fly)'-(TTL_timestamps.fly(1)-Timestamps_of_first_samples_usec(1))/1e6;    RPL_fly.a = RPL_a(idx_fly)';    RPL_fly.d = RPL_d(idx_fly)';    RPL_fly.r = RPL_r(idx_fly)';
    RPL_pst.t = RPL_t(idx_pst)'-(TTL_timestamps.pst(1)-Timestamps_of_first_samples_usec(1))/1e6;    RPL_pst.a = RPL_a(idx_pst)';    RPL_pst.d = RPL_d(idx_pst)';    RPL_pst.r = RPL_r(idx_pst)';
    save([analysis_directory,'\Ripples',CSC_file.name(6:end-4),'_AF.mat'],'RPL_pre','RPL_fly','RPL_pst','P_ave');