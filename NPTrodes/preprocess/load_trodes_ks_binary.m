% function D=load_trodes_ks_binary(jsonFile, type, index, varargin)

Fs = 30000;
Fs_lfp = 1500;

%% load kilosort data
% ks_file = dir(fullfile(pwd,'*kilosort','*.probe*'));
ts_file = dir(fullfile(pwd,'*kilosort','*timestamps.dat'));
% fname_1 = fullfile(ks_file(1).folder,ks_file(1).name);
% fname_2 = fullfile(ks_file(2).folder,ks_file(2).name);
tsname = fullfile(ts_file.folder,ts_file.name);
temp_file = dir(fullfile(pwd,'*','temp_wh_probe*'));
tempname_1 = fullfile(temp_file(1).folder,temp_file(1).name);
tempname_2 = fullfile(temp_file(2).folder,temp_file(2).name);

% fname = ops.fname_1;
num_chan = 200; % number of channels
% fname = fname_1;
fname = tempname_1;
file = dir(fname);
samples=file.bytes/2/num_chan;
samples/30000/60/60

%% load temp_wh.dat
D.Data=memmapfile(fname,'Format',{'int16' [num_chan samples] 'mapped'});
ks_data = double(D.Data.Data(1).mapped(1,1:30000*60*10));

% Filter temp_wh.dat
lfp_ks = bandpass(ks_data,[0.5,300],30000,ImpulseResponse="iir")';
% [b,a] = butter(3,[0.5/Fs,300/Fs],'bandpass');
% lfp_ks = filter(b,a,ks_data)';
ap_ks = highpass(ks_data,600,30000)';

% load timestamps
% ts_ks = readTrodesExtractedDataFile(tsname);
% ts_ks = double(ts_ks.fields.data);
% samples_true = length(ts_ks);
% samples_diff = samples-samples_true;

ts_ap = 1/Fs:1/Fs:60*10; % timstamps for 30kHz
ts_lfp = 1/Fs_lfp:1/Fs_lfp:60*10; % timestamps for 1500Hz

% dowsample to 1500Hz
lfp_ks_1500 = interp1(ts_ap,lfp_ks,ts_lfp)';   % 0.85 steepness, 60 dB

%% load spikes extracted by spikegadgets
ap_file = dir(fullfile(pwd,'*spikeband','*nt1*'));
ap_data = readTrodesExtractedDataFile(fullfile(ap_file(1).folder,ap_file(1).name));
ap = double(ap_data.fields.data(1:30000*60*10));
ap = highpass(ap,600,30000);

%% Load lfp extracted by Trodes
Fs_lfp = 1500;
lfp_file = dir(fullfile(pwd,'*LFP','*nt1*'));
lfp_data = readTrodesExtractedDataFile(fullfile(lfp_file(1).folder,lfp_file(1).name));
lfp = double(lfp_data.fields.data(1:Fs_lfp*60*10));
lfp = bandpass(lfp,[0.5 300],Fs_lfp);
% [b,a] = butter(3,[0.5/Fs_lfp,300/Fs_lfp],'bandpass');
% lfp2 = filter(b,a,lfp);

%% Lags b/w kilosort and spikeband
% [r,lags] = xcorr(ap(10:end),ap_ks(1:length(ap_ks)-9),50000);
[r,lags] = xcorr(ap,ap_ks,50000);
lag_est = lags(r == max(r));
ap_corrected = ap(1+lag_est:end);
ap_ks_corrected = ap_ks(1:length(ap_ks)-lag_est);

[r_corrected,lags_corrected] = xcorr(ap_corrected,ap_ks_corrected,50000);

figure; plot(lags,r,'b','DisplayName','temp-wh vs spikeband'); hold on
plot(lags_corrected,r_corrected,'r','DisplayName','Corrected')
plot([0,0],[min(r)*1.5,max(r)*1.5],'k','DisplayName','lag = 0')
hold off
xlim([-100,100])
ylim([min(r)*1.5,max(r)*1.5])
legend()

%% Lags b/w kilosort and lfp
% [r2,lags2] = xcorr(lfp,lfp_ks_1500,1000);
% lag_est2 = lags2(r2 == max(r2));
% lfp_corrected = lfp(1+lag_est2:end);
% lfp_ks_corrected = lfp_ks_1500(1:length(lfp_ks_1500)-lag_est2);
% 
% [r_corrected,lags_corrected] = xcorr(lfp_corrected,lfp_ks_corrected,1000);
% 
% figure; plot(lags2,r2,'b','DisplayName','temp-wh vs spikeband'); hold on
% plot(lags_corrected,r_corrected,'r','DisplayName','Corrected')
% plot([0,0],[min(r2)*1.5,max(r2)*1.5],'k','DisplayName','lag = 0')
% hold off
% xlim([-100,100])
% ylim([min(r2)*1.5,max(r2)*1.5])
% legend()

%% Compare temp_wh.dat vs spikeband
%=== Show a sample of the session
t1 = 50;   t2 = t1+0.01;
figure('units','normalized','outerposition',[0 0 1 1]);
tl = tiledlayout(5,2,'TileSpacing','tight');
for i=1:10
    nexttile;
    % plot(ts_ks(ts_ks>i*t1 & ts_ks<i*t2),ks_data(ts_ks>i*t1 & ts_ks<i*t2),'k');   hold on;
    plot(ts_ap(ts_ap>i*t1 & ts_ap<i*t1+0.01),ap_ks_corrected(ts_ap>i*t1 & ts_ap<i*t1+0.01),'k'); hold on
    plot(ts_ap(ts_ap>i*t1 & ts_ap<i*t1+0.01),ap_corrected(ts_ap>i*t1 & ts_ap<i*t1+0.01),'r');
    % ylim([-1e3 1e3]);
    xlabel('Time (s)'); ylabel('Amplitude');
    hold off
end

title(tl,'temp-wh.dat vs spikeband extracted by Trodes: red = Trodes, black = temp-wh.dat')


%% Compare temp_wh.dat vs lfp band
t1 = 50;   t2 = 10;
figure('units','normalized','outerposition',[0 0 1 1]);
tl = tiledlayout(5,2,'TileSpacing','tight');
for i=1:10
    nexttile;
    % plot(ts_ks(ts_ks>i*t1 & ts_ks<i*t2),ks_data(ts_ks>i*t1 & ts_ks<i*t2),'k');   hold on;
    plot(ts_lfp(ts_lfp>i*t1 & ts_lfp<i*t1+t2),lfp_ks_1500(ts_lfp>i*t1 & ts_lfp<i*t1+t2),'k'); hold on
    plot(ts_lfp(ts_lfp>i*t1 & ts_lfp<i*t1+t2),lfp(ts_lfp>i*t1 & ts_lfp<i*t1+t2),'r');
    % ylim([-1e3 1e3]);
    xlabel('Time (s)'); ylabel('Amplitude');
    hold off
end

title(tl,'temp-wh.dat vs LFP extracted by Trodes: red = Trodes, black = temp-wh.dat')


%% load raw.dat file
raw_file = dir(fullfile(pwd,'*raw','*raw*'));
raw_probe1 = fullfile(raw_file(1).folder,raw_file(1).name);
raw_probe2 = fullfile(raw_file(2).folder,raw_file(2).name);

%%
clc
fid = fopen(raw_probe1,'rb');
% header = fread(fid,[num_chan 150000],"uint8");
data_format = 'uint8';
header = fread(fid,1000,data_format);
header = char(header')
fieldsLoc  = strfind(header,'First_timestamp');
fseek(fid, fieldsLoc+25, -1);
values = fread(fid,[200,150000],[data_format, '=>int16']);
% frewind(fid)
% figure; plot(header(1,:))

%%
num_chan = 200; % number of channels
fname = raw_probe1;
file = dir(fname);
samples=file.bytes/2/num_chan;
samples/30000/60/60

D2.Data=memmapfile(fname,'Format',{'int16' [num_chan samples] 'mapped'});
raw_data = double(D.Data.Data(1).mapped(1,1:30000*60*10));

% %% Spectrogram
% nsc = 120; % window size is 120 sample (120 ms)
% nov = nsc/2; % 50% overlap
% Fs = 30000;
% % [mag_spctr,f_spctr,t_spectr] = spectrogram(LFP,hamming(nsc),nov,Fs);
% % [mag_spctr,f_spctr,t_spctr] = spectrogram(LFP,nsc,nov,1024,Fs,'centered');
% [mag_spctr,f_spctr,t_spctr] = spectrogram(lfp,nsc,nov,1024,Fs);
% % [mag_spctr,f_spctr,t_spectr] = spectrogram(LFP,hamming(nsc),nov);
% % f_spctr = f_spctr*Fs/2; % frequency (Hz)
% figure('units','normalized','outerposition',[0 0.3 1 0.5]);
% imagesc(t_spctr*2/60,f_spctr(f_spctr<450),abs(mag_spctr(f_spctr<450,:)))
% colorbar
% % title(sprintf('LFP (1-150Hz) spectrogram batdate = %s',batdate))
% set(gca,'YDir','normal')
% ylabel('Frequency (Hz)')
% xlabel("Time (min)")
% 
% %% normalize spectrogram
% % normalize spectrogram
% S = abs(mag_spctr(f_spctr<150,:));
% f_spctr = f_spctr(f_spctr<150);
% norm_mad = mad(S,1,2);
% norm_med = median(S,2);
% S_norm = (S -norm_med) ./ norm_mad;
% M = mean(S_norm,1);
% % visualize normalized spectrogram
% figure('units','normalized','outerposition',[0 0.3 1 0.5]);
% imagesc(t_spctr*2/60,f_spctr,S_norm)
% colorbar
% % title(sprintf('LFP (1-150Hz) normalized spectrogram batdate = %s',batdate))
% set(gca,'YDir','normal')
% ylabel('Frequency (Hz)')
% xlabel("Time (min)")
% 
% 
% 
% 
% %%
% % headerText = fread(fid,10000,'uint8');
% % headerText = char(headerText')
% headerText = fread(fid,100,'uint8');
% headerText = double(headerText');