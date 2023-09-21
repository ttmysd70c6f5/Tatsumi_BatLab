% load behavior data
beh_dir = dir(fullfile(pwd,'Ext_Behavior*','Extracted_Behavior*'));
load(fullfile(beh_dir.folder,beh_dir.name),'a','a_abs','angle','bat_clr','bat_nms','bat_pair_nms','bat_pairs',...
    'batdate','bflying','Fs','f_num','Group_name','n_tags','r','r_lim','t','T','v','v_abs','v_th','wBeats');        % Processed Behavioral Data
% load ephys data
ephys_dir = dir(fullfile(cd,'Ext_Ephys*','CSC_Analysis*','LFP_Ext*'));
ephys_data = struct([]);
for bat = 1:length(ephys_dir)
    load(fullfile(ephys_dir(bat).folder,ephys_dir(bat).name),'LFP_data')
    ephys_data(bat).LFP = LFP_data;
end

%% Figure settings

%% Extract when implanted bats are close
% Extract when both implanted bats are resting
nimpbat = length(ephys_data);
imp_bat_pairs = nchoosek(1:nimpbat,2);
imp_bats = false(1,n_tags);
for bat = 1:nimpbat
    imp_bats(ephys_data(bat).LFP.imp_bat) = true;
end

imp_resting = ~any(bflying(:,imp_bats),2);
r_imp_resting = r(imp_resting,:,imp_bats);
t_imp_resting = t(imp_resting);
d_imp_resting = zeros([size(r_imp_resting,[1,2]),size(imp_bat_pairs,1)]);
for i = size(imp_bat_pairs,1)
    d_imp_resting = abs(r_imp_resting(:,:,imp_bat_pairs(1)) - r_imp_resting(:,:,imp_bat_pairs(2)));
end
d_imp_resting = vecnorm(d_imp_resting,2,2);

% plot the distance between implanted bats while resting
figure('units','normalized','outerposition',[0 0.3 1 0.5]);
plot(t_imp_resting,d_imp_resting)

%% Spectrogram
bat = 1;
TT1 = 1;
Freq_band = [0.5 150];
Fs = 1000;
% LFP = ephys_data(bat).LFP.LFP(TT1,1e5+1:end-1e5);
% t_1KHz = ephys_data(bat).LFP.t_1KHz(TT1,1e5+1:end-1e5);
LFP = ephys_data(bat).LFP.LFP(TT1,:);
t_1KHz = ephys_data(bat).LFP.t_1KHz(TT1,:);

LFP = bandpass(LFP,Freq_band,Fs);

% nsc = floor(length(LFP)/120); % window size is 120 sample (120 ms)
% nov = floor(nsc/2); % 50% overlap
nsc = 120; % window size is 120 sample (120 ms)
nov = nsc/2; % 50% overlap
% [mag_spctr,f_spctr,t_spectr] = spectrogram(LFP,hamming(nsc),nov,Fs);
% [mag_spctr,f_spctr,t_spctr] = spectrogram(LFP,nsc,nov,1024,Fs,'centered');
[mag_spctr,f_spctr,t_spctr] = spectrogram(LFP,nsc,nov,1024,Fs);
% [mag_spctr,f_spctr,t_spectr] = spectrogram(LFP,hamming(nsc),nov);
% f_spctr = f_spctr*Fs/2; % frequency (Hz)
figure('units','normalized','outerposition',[0 0.3 1 0.5]);
imagesc(t_spctr*2/60,f_spctr(f_spctr<150),log(abs(mag_spctr(f_spctr<150,:))))
colorbar
title(sprintf('LFP (1-150Hz) spectrogram batdate = %s',batdate))
set(gca,'YDir','normal')
ylabel('Frequency (Hz)')
xlabel("Time (min)")

%% normalize spectrogram
% normalize spectrogram
S = abs(mag_spctr(f_spctr<150,:));
f_spctr = f_spctr(f_spctr<150);
norm_mad = mad(S,1,2);
norm_med = median(S,2);
S_norm = (S -norm_med) ./ norm_mad;
M = mean(S_norm,1);
% visualize normalized spectrogram
figure('units','normalized','outerposition',[0 0.3 1 0.5]);
imagesc(t_spctr*2/60,f_spctr,S_norm)
colorbar
title(sprintf('LFP (1-150Hz) normalized spectrogram batdate = %s',batdate))
set(gca,'YDir','normal')
ylabel('Frequency (Hz)')
xlabel("Time (min)")

%% PCA for spectrogram
[coeff_pca,score,latent,~,explained,mu] = pca(S_norm');
var_combined = sum(latent(1:2));
fs_range = 1:149;
var_div = zeros(1,length(fs_range));
for i = 1:length(fs_range)
    div_fs = fs_range(i);
    var_div(i) = sum(var(S_norm(f_spctr<=div_fs,:),0,2)) + sum(var(S_norm(f_spctr>div_fs,:),0,2));
    var_div(i) = var_div(i) / var_combined;
end

% proportion of variance
sz = 40;
figure;
hold on
scatter(1,explained(1),sz,'r','filled')
scatter(2,explained(2),sz,'b','filled')
scatter(3:length(explained),explained(3:end),sz,'k','filled')
hold off
xlabel('Principal component')
ylabel('Propotion of variance explained')

%% PCA score
figure;
scatter(score(:,1),score(:,2))

%% Relative normalized power
figure; hold on
for i = 1:100
    plot(f_spctr,coeff_pca(:,i),'DisplayName',sprintf('PC%d',i));
end
hold off
legend()
ylabel('Relative normalized power')
xlabel('Frequency (Hz)')

%% Two freq band normalized LFP power
freq_band_pc = [1 75; 75 150];

figure('units','normalized','outerposition',[0 0.2 1 0.7]);
tiledlayout(2,1,'TileSpacing','tight')
for i = 1:2
    nexttile; hold on
    for bat = 1:nimpbat
        LFP = ephys_data(bat).LFP.LFP(1,1e5+1:end-1e5);
        LFP = bandpass(LFP,freq_band_pc(i,:),1000);
        t_1KHz = ephys_data(bat).LFP.t_1KHz(1e5+1:end-1e5);
        plot(t_1KHz/60,LFP)
        title(sprintf('LFP traces frequency band = [%d %d] batdate = %s',freq_band_pc(i,1),freq_band_pc(i,2),batdate))
        ylabel('Amplitude')
        ylim([-4000 4000])
    end
    hold off
end
xlabel('Time (min)')

%% LFP correlation for each frequency band
corr_lfp_two_band = zeros(1,2);
for i = 1:2
    LFP_1 = ephys_data(1).LFP.LFP(1,1e5+1:end-1e5);
    LFP_1 = bandpass(LFP_1,freq_band_pc(i,:),1000);
    LFP_2 = ephys_data(2).LFP.LFP(1,1e5+1:end-1e5);
    LFP_2 = bandpass(LFP_2,freq_band_pc(i,:),1000);
    c = corrcoef(LFP_1,LFP_2);
    corr_lfp_two_band(i) = c(1,2);
end

%%
LFP_1 = ephys_data(1).LFP.LFP(1,1e5+1:end-1e5);
LFP_1 = bandpass(LFP_1,Freq_band,1000);
LFP_2 = ephys_data(2).LFP.LFP(1,1e5+1:end-1e5);
LFP_2 = bandpass(LFP_2,Freq_band,1000);
c = corrcoef(LFP_1,LFP_2);
c = c(1,2);

%% Preprocessing lfp data
% High (30-150Hz)and low (1-29Hz) frequency
% bandstop filter
for bat = 1:nimpbat
    for TT = 1:size(ephys_data(bat).LFP.LFP,1)
        raw_voltage = ephys_data(bat).LFP.LFP(TT,:);
        filtered_voltage = bandstop(raw_voltage,[59.5,60.5],1000);
        filtered_voltage = bandstop(filtered_voltage,[119.5,120.5],1000);
        ephys_data(bat).LFP.LFP_filtered(TT,:) = filtered_voltage;
    end
end

%=== Show a sample of the session
for bat = 1:nimpbat
    t_1KHz = ephys_data(bat).LFP.t_1KHz;
    LFP = ephys_data(bat).LFP.LFP;
    LFP_filtered = ephys_data(bat).LFP.LFP_filtered;
    for TT = 1:size(LFP,1)
        t1 = 150;   t2 = t1+10;
        figure('units','normalized','outerposition',[0 0 1 1]);
        tiledlayout(5,2,'TileSpacing','tight');
        for i=1:10
            nexttile;
            plot(t_1KHz(t_1KHz>i*t1 & t_1KHz<i*t1+10),LFP(TT,t_1KHz>i*t1 & t_1KHz<i*t1+10),'k');
            hold on;
            plot(t_1KHz(t_1KHz>i*t1 & t_1KHz<i*t1+10),LFP_filtered(TT,t_1KHz>i*t1 & t_1KHz<i*t1+10),'r');
            % ylim([-1e3 1e3]);
            xlabel('Time (s)'); ylabel('Voltage (uV)');
            hold off
        end
    end
end
sgtitle(['LFP, TT', num2str(TT), ' , ', num2str(batdate)]);

%% LFP traces of two implanted bats
% plot the distance between implanted bats while resting
figure('units','normalized','outerposition',[0 0.3 1 0.5]);
hold on
for bat = 1:nimpbat
    LFP = ephys_data(bat).LFP.LFP(1,:);
    t_1KHz = ephys_data(bat).LFP.t_1KHz;
    plot(t_1KHz,LFP)
end
hold off
ylim([-4000 4000])

%% Moving pearson correlation
freq_band_pc = [1 50; 100 150];
for freq = 1:2
    wd = 1*60; % window size (sec). must be < 2*1e6
    lag = 5; % time lag (sec)
    TT1 = 4;
    TT2 = 4;
    LFP_imp1 = ephys_data(1).LFP.LFP(TT1,:);
    LFP_imp1 = bandpass(LFP_imp1,freq_band_pc(freq,:),1000);
    LFP_imp2 = ephys_data(2).LFP.LFP(TT2,:);
    LFP_imp2 = bandpass(LFP_imp2,freq_band_pc(freq,:),1000);
    t_1KHz_imp1 = ephys_data(1).LFP.t_1KHz;
    t_1KHz_imp2 = ephys_data(2).LFP.t_1KHz;
    tstamp_corr = 1:lag*1e3:length(t_1KHz_imp1);
    t_corr = t_1KHz_imp1(tstamp_corr);
    lfp_corr = zeros(1,length(tstamp_corr));
    for i = 1:length(tstamp_corr)
        ts = tstamp_corr(i);
        l = wd/2*1e3;
        c = corrcoef(LFP_imp1(max(ts-l,1):min(ts+l,length(LFP_imp1))),LFP_imp2(max(ts-l,1):min(ts+l,length(LFP_imp2))));
        lfp_corr(i) = c(1,2);
    end
    
    figure('units','normalized','outerposition',[0 0.2 1 0.7]);
    tiledlayout(3,1,'TileSpacing','tight')
    ax(1) = nexttile; plot(t_corr/60,lfp_corr);
    title(sprintf('LFP corr: Freq band = [%d %d], Date = %s, window = %dmin, lag = %dmin',freq_band_pc(freq,1),freq_band_pc(freq,2), batdate,wd/60,lag/60))
    ylim([-0.3 0.3])
    ylabel('Correlation')
    ax(2) = nexttile;
    % plot(t_1KHz_imp1(1e6:end)/60,LFP_imp1(1,1e6:end)); hold on
    % plot(t_1KHz_imp2(1e6:end)/60,LFP_imp2(1,1e6:end))
    plot(t_1KHz_imp1/60,LFP_imp1); hold on
    plot(t_1KHz_imp2/60,LFP_imp2)
    hold off
    ylabel('LFP traces')
    ylim([-4000 4000])
    
    ax(3) = nexttile;
    plot(t_imp_resting/60,d_imp_resting)
    xlabel('Time (min)')
    ylabel('distance')
    linkaxes(ax,'x')
end