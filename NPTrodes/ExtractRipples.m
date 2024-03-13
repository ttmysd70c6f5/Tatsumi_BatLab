function Ripples = ExtractRipples(lfp,fig_output)
%% Experiment info
gain = 50; % LFP gain
fs = 1500; % LFP sampling rate

%% Load data
% lfp_file = dir(fullfile(lfp_dir,'*LFP_nt*.dat')); % list of LFP .dat files
% lfp_data = readTrodesExtractedDataFile(fullfile(lfp_file(ch).folder, lfp_file(ch).name)); % Load LFP data
% lfp = double(lfp_data.fields.data);
% % str = extract(lfp_file(ch).name,digitsPattern);
% % ch_depth = 10*(str2double(str{3})-1000);
% 
% % Subtract median
% load(fullfile(lfp_file(1).folder,'lfpmedian.mat'),'lfpmed1')
% lfp = lfp - lfpmed1;
% lfp = lfp - median(lfp);
% 
% lfp = lowpass(lfp,400,fs);
% 
% % lfp = lfp / gain; % convert to µV

%% Putative ripples
% Extract ripple (80-160Hz)
lfp_ripple = bandpass(lfp',[80 160],fs)'; % Extract high-frequency ripples
ripple_power = abs(hilbert(lfp_ripple')'); % Hilbert transform of lfp
ripple_pwmean = mean(ripple_power,2);
ripple_pwsd = std(ripple_power,0,2); % std of lfp power

th = 7;
Ripples = ripple_power > ripple_pwmean + th*ripple_pwsd;
Ripples = find(diff([Ripples(1), Ripples]) == 1);


%% Fig

wd = 0.1; % window size in ms
ripple_wv = zeros(length(Ripples),fs*2*wd+1);

idx_pre = 0; % initialize parameter
ripple_cnt = 0; % # of ripples
to_del = false(size(Ripples));
for i=1:length(Ripples)
    idx_now = Ripples(i);
    
    if idx_now > fs*wd && idx_now - idx_pre > fs*wd && idx_now < length(lfp) - fs*wd
%     if idx_now > fs*wd && idx_now < length(lfp) - fs*wd
%         nexttile
%         plot(-fs*0.1:fs*0.1, lfp_ripple(ch,Ripples_idx(i)-fs*0.1:Ripples_idx(i)+fs*0.1),'k')
%         title([num2str(Ripples_idx(i)/1500),' (sec)'])
    
        ripple_wv(i,:) = lfp_ripple(Ripples(i)-fs*wd:Ripples(i)+fs*wd);
        ripple_cnt = ripple_cnt + 1;
    else
        to_del(i) = true;
    end

    idx_pre = idx_now;
end

ripple_wv(to_del,:) = [];
Ripples(to_del) = [];


if fig_output
    % figure
    plot([-wd:1/fs:wd].*1000,mean(ripple_wv),'k')
%     title(sprintf('putative ripple: electrode %d',ch_depth/10))
    title('putative ripple')
    ylabel('Amplitude')
    xlabel('Time (ms)')
    % xlim([-wd*fs wd*fs])
end