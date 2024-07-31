% function InterBrainCoupling(recID_1,recID_2)
% Inter-brain coupling analysis
%% General exp info
% recName = mtData.recName;

%% Load neural data
recListPath = 'H:\My Drive\YartsevLab\Experiment\Multi_reward\Log\Multireward_exp_list_2024.xlsx'; % table of rec info
allRecs = GetRecList(recListPath); % load the rec info list

% Load data
% recID_1 = 3; recID_2 = 4; % 240207, #14633 and 14462, MEC
% recID_1 = 6; recID_2 = 7; % 240208, #14633 and 14462, MEC
% recID_1 = 3; recID_2 = 4; % 240209, #14633 and 14462, MEC
recID_1 = 12; recID_2 = 13; % 240210_1, #14633 and 14462, MEC
% recID_1 = 3; recID_2 = 4; % 240212, #14633 and 14462, MEC
dataTypes = 'barslm'; % data types to be loaded
[BhvData_1,AnalogData_1,RewardData_1,spikeData_1,lfpData_1,mtData_1] = LoadExperimentData(recID_1,dataTypes);
[BhvData_2,AnalogData_2,RewardData_2,spikeData_2,lfpData_2,mtData_2] = LoadExperimentData(recID_2,dataTypes);


disp('Extracting flights ...')
[Flights_1, ~] = extractFlights(BhvData_1,RewardData_1,mtData_1,'include','all');
[Flights_2, ~] = extractFlights(BhvData_2,RewardData_2,mtData_2,'include','all');

% Ehys bats
bephys_1 = 6;
bephys_2 = 7;

%%
% Channel in the MEC (Caution! The depth info in the loaded LFP structure
% is incorrect. A smaller channel number corresponds to a deeper region)
recName = mtData_1.recName;
switch recName
    case '240207'
        chan_1 = 46;
        chan_2 = 9;
    case '240208'
        chan_1 = 46; % MEC for bat 1: tip channel
        chan_2 = 7; % MEC for bat 2: tip channel
    case '240210_1'
        chan_1 = 49;
        chan_2 = 9;
end
pr_1 = 1; % probe number
pr_2 = 1; % probe number

lfp_raw_1 = double(lfpData_1.lfp{pr_1}(chan_1,:));
t_lfp_1 = lfpData_1.global_sample_timestamps_usec/1e6;

lfp_raw_2 = double(lfpData_2.lfp{pr_2}(chan_2,:));
t_lfp_2 = lfpData_2.global_sample_timestamps_usec/1e6;

%% Align two recordings
% t_bhv_max = max([BhvData_1.t(end), BhvData_2.t(end)])
t_bhv_idx_max = min([length(BhvData_1.t), length(BhvData_2.t)]);
t_bhv_max = min([BhvData_1.t(end), BhvData_2.t(end)]);
t_ephys_max = min([max(spikeData_2.allGlobalSpikeTimes_usec), max(spikeData_1.allGlobalSpikeTimes_usec)] / 1e6);

t_min = 0;
t_max = min([t_bhv_max,t_ephys_max]);

sr_lfp = 2500; % sampling rate
t_align = t_min:1/sr_lfp:t_max; % global time stamp to align

% idx = t_lfp_1 > 0 & t_lfp_1 < t_max;
% lfp_1 = lfp_1(idx); t_lfp_1 = t_lfp_1(idx);
[~,idx_1,~] = unique(t_lfp_1);
lfp_align_1 = interp1(t_lfp_1(idx_1),lfp_raw_1(idx_1),t_align,'linear');

% idx = t_lfp_2 > 0 & t_lfp_2 < t_max;
% lfp_2 = lfp_2(idx); t_lfp_2 = t_lfp_2(idx);
[~,idx_2,~] = unique(t_lfp_2);
lfp_align_2 = interp1(t_lfp_2(idx_2),lfp_raw_2(idx_2),t_align,'linear');

% Visualize the alignment
fig = figure;
set(fig,'units','normalize','Position',[0 0 1 1])
subplot(2,1,1)
hold on
plot(t_lfp_1,lfp_raw_1,'-r')
plot(t_align,lfp_align_1,'-k')
hold off
xlim([t_min t_max])
title('LFP: Bat 1')

subplot(2,1,2)
hold on
plot(t_lfp_2,lfp_raw_2,'-b')
plot(t_align,lfp_align_2,'-k')
hold off
xlim([t_min t_max])
title('LFP: Bat 2')

lfp_1 = lfp_align_1; lfp_2 = lfp_align_2; t_lfp = t_align;

%% Define behavioral state
n_bats = BhvData_1.n_tags;
T = BhvData_1.T;
r = BhvData_1.r;
bflying = BhvData_1.bflying;
blanding = 1 - bflying;

% params
th_social = 0.5;
th_nbats = 2;

% Inter-bat distance
d = NaN(T,n_bats,n_bats);
for bb = 1:n_bats
    for bbb = 1:n_bats
        d(:,bb,bbb) = vecnorm(r(:,:,bb)-r(:,:,bbb),2,2); % [T,B,B]
    end
end

% Count # of neighor bats
n_neighbors = zeros(T,n_bats); % [T,B]
for bb = 1:n_bats
    n_neighbors(:,bb) = sum(d(:,[1:bb-1,bb+1:end],bb)<th_social,2); % [T,1]
end

% Define bat states
% 0 = flying, -1 = individual resting (no neighor bat), -2 = resting with the target bat, 1-9 = # of neighbors
bstate = NaN(T,n_bats);
for bb = 1:n_bats
    % flying
    bstate(bflying(:,bb)==1,bb) = 0;
    % Resting with non-ephys bat
    bstate(blanding(:,bb)==1 & n_neighbors(:,bb)>0 & n_neighbors(:,bb)<th_nbats,bb) = 1; % # of neighbors
    % Resting with non-ephys bat
    bstate(blanding(:,bb)==1 & n_neighbors(:,bb)>=th_nbats,bb) = 1; % # of neighbors
    % Alone resting
    bstate(blanding(:,bb)==1 & n_neighbors(:,bb)==0,bb) = -1;
    % Resting with other ephys bat
    if bb == bephys_1
        bstate(blanding(:,bb)==1 & d(:,bb,bephys_2)<th_social,bb) = -2;
    elseif bb == bephys_2
        bstate(blanding(:,bb)==1 & d(:,bb,bephys_1)<th_social,bb) = -2;
    end
end



%% Align behavioral state to LFP recordings
t = BhvData_1.t;
n_bats = BhvData_1.n_tags;
bstate_new = NaN(length(t_lfp),n_bats);

for bb = 1:n_bats
    bstate_new(:,bb) = interp1(t,bstate(:,bb),t_lfp,'nearest');
end

bstate = bstate_new;

%% Compute raw LFP spectrogram among session
Freq_band = [0.5 200];
sr_lfp = 2500; % sampling rate
th_lfp = 150;

lfp_filt_1 = bandpass(lfp_1,Freq_band,sr_lfp);
lfp_filt_2 = bandpass(lfp_2,Freq_band,sr_lfp);

nsc = 300; % window size (120 ms)
nov = nsc/2; % 50% overlap
% [mag_spctr,f_spctr,t_spectr] = spectrogram(LFP,hamming(nsc),nov,Fs);
% [mag_spctr,f_spctr,t_spctr] = spectrogram(LFP,nsc,nov,1024,Fs,'centered');
[mag_spctr_1,f_spctr_1,t_spctr_1] = spectrogram(lfp_filt_1,nsc,nov,1024,sr_lfp); S_abs_1 = abs(mag_spctr_1(f_spctr_1<th_lfp,:)); f_spctr_1 = f_spctr_1(f_spctr_1<th_lfp,:);
[mag_spctr_2,f_spctr_2,t_spctr_2] = spectrogram(lfp_filt_2,nsc,nov,1024,sr_lfp); S_abs_2 = abs(mag_spctr_2(f_spctr_2<th_lfp,:)); f_spctr_2 = f_spctr_2(f_spctr_2<th_lfp,:);
% [mag_spctr,f_spctr,t_spectr] = spectrogram(LFP,hamming(nsc),nov);
% f_spctr = f_spctr*Fs/2; % frequency (Hz)

%% Normalization of LFP spectrogram
% peak normalization
% S_peak_1 = max(S_abs_1,[],2); S_peak_2 = max(S_abs_2,[],2);
% S_norm_1 = S_abs_1 ./ S_peak_1; S_norm_2 = S_abs_2 ./ S_peak_2;

% z-score normalization
S_norm_1 = normalize(S_abs_1,2,'zscore'); S_norm_2 = normalize(S_abs_2,2,'zscore');

%% PCA for spectrogram
[coeff_pca_1,score_1,latent_1,~,explained_1,mu_1] = pca(S_norm_1');
[coeff_pca_2,score_2,latent_2,~,explained_2,mu_2] = pca(S_norm_2');

%% LFP correlation
state_list = [-2,-1,0,1];
corr_lfp = NaN(2,length(state_list));

sr_lfp = 2500;
Freq_band_b1 = [0.5,25];
Freq_band_b2 = [25 150];
lfp_b1_1 = bandpass(lfp_1,Freq_band_b1,sr_lfp); lfp_b1_2 = bandpass(lfp_2,Freq_band_b1,sr_lfp);
lfp_b2_1 = bandpass(lfp_1,Freq_band_b2,sr_lfp); lfp_b2_2 = bandpass(lfp_2,Freq_band_b2,sr_lfp);

x = bstate(:,6);
for i = 1:length(state_list)
    corr_now = corrcoef(lfp_b1_1(x==state_list(i)), lfp_b1_2(x==state_list(i)));
    corr_lfp(1,i) = corr_now(1,2);
    corr_now = corrcoef(lfp_b2_1(x==state_list(i)), lfp_b2_2(x==state_list(i)));
    corr_lfp(2,i) = corr_now(1,2);
end

% Sliding correlation
lfp_movcorr_1 = movcorr(lfp_b1_1',lfp_b1_2',sr_lfp*30);
lfp_movcorr_2 = movcorr(lfp_b2_1',lfp_b2_2',sr_lfp*30);
lfp_movcorr = {lfp_movcorr_1,lfp_movcorr_2};

corr_lfp

% Save lfp corr
recName = mtData_1.recName;
figDir = 'F:\Data\GroupForaging\Analysis\figure';
saveDir = fullfile(figDir,'InterBrainCoupling');
if isempty(dir(saveDir)); mkdir(saveDir); end
save(fullfile(saveDir,sprintf('lfpCorr_%s.mat',recName)),'corr_lfp','lfp_movcorr')


%% Summary across sessions
figDir = 'F:\Data\GroupForaging\Analysis\figure';
saveDir = fullfile(figDir,'InterBrainCoupling');

Files = dir(fullfile(saveDir,'lfpCorr*.mat'));
lfpCorrs = cell(length(Files),1);
for i = 1:length(Files)
    load(fullfile(Files(i).folder,Files(i).name),'corr_lfp')
    lfpCorrs{i} = corr_lfp;
end
%% FIGURES
%%%% Visualize behavioral states

x = bstate(:,7);
[x_start,x_end,x_list] = sectionalize(x);
c = [0.6 0.2 0.2; % red: resting with the target bat
    0.3 0.7 0.9; % light blue: individual resting
    0 0 0; % flying
    0.1 0.6 0.2]; % green: resting with non-target bat
    
figure
set(gcf,'Units','Normalized','Position',[0 0.3 1 0.4])
set(gca,'Units','Normalized','Position',[0.05 0.1 0.9 0.8])
hold on
for i = 1:length(x_list)
% for i = 4
    xregion(x_start{i},x_end{i},'FaceColor',c(i,:),'FaceAlpha',0.5)
end
hold off
xlim([0 T])
% xlim([t_vector(1),t_vector(end)])


%%%% Visualize spectrogram
%%% FIGURE: bat 1
figure
set(gcf,'Units','Normalized','Position',[0 0.3 1 0.5])
set(gca,'Units','Normalized','Position',[0 0.1 1 0.8])
% Behavioral state
x = bstate(:,6);
[x_start,x_end,x_list] = sectionalize(x);
c = [0.6 0.2 0.2; % red: resting with the target bat
    0.3 0.7 0.9; % light blue: individual resting
    0 0 0; % flying
    0.1 0.6 0.2]; % green: resting with non-target bat

ax(1) = subplot(3,1,1);
hold on
for i = 1:length(x_list)
% for i = 4
    xregion(t_lfp(x_start{i}),t_lfp(x_end{i}),'FaceColor',c(i,:),'FaceAlpha',0.5)
end
hold off
xlim([t_min t_max])
yticks([])
% LFP spectrogram
ax(2) = subplot(3,1,2);
imagesc(t_spctr_1,f_spctr_1,log(S_abs_1))
clb = colorbar; clb.Location = 'east';
title('Bat1: LFP spectrogram')
set(gca,'YDir','normal')
ylabel('Frequency (Hz)')
xlabel("Time (sec)")
xlim([t_min,t_max])

ax(3) = subplot(3,1,3);
imagesc(t_spctr_1,f_spctr_1,S_norm_1)
clb = colorbar; clb.Location = 'east';
title('Bat1: Normalized spectrogram')
set(gca,'YDir','normal')
ylabel('Frequency (Hz)')
xlabel("Time (sec)")
xlim([t_min,t_max])


linkaxes(ax,'x')

%%% FIGURE: bat 2
figure('units','normalized','outerposition',[0 0.3 1 0.5]);

% Behavioral state
x = bstate(:,7);
[x_start,x_end,x_list] = sectionalize(x);
c = [0.6 0.2 0.2; % red: resting with the target bat
    0.3 0.7 0.9; % light blue: individual resting
    0 0 0; % flying
    0.1 0.6 0.2]; % green: resting with non-target bat

ax(1) = subplot(3,1,1);
hold on
for i = 1:length(x_list)
% for i = 4
    xregion(t_lfp(x_start{i}),t_lfp(x_end{i}),'FaceColor',c(i,:),'FaceAlpha',0.5)
end
hold off
xlim([t_min t_max])
yticks([])
% LFP spectrogram
ax(2) = subplot(3,1,2);
imagesc(t_spctr_2,f_spctr_2,log(S_abs_2))
clb = colorbar; clb.Location = 'east';
title('Bat2: LFP spectrogram')
set(gca,'YDir','normal')
ylabel('Frequency (Hz)')
xlabel("Time (sec)")
xlim([t_min,t_max])
% Normalized LFP spectrogram
ax(3) = subplot(3,1,3);
imagesc(t_spctr_2,f_spctr_2,S_norm_2)
clb = colorbar; clb.Location = 'east';
title('Bat2: Normalized spectrogram')
set(gca,'YDir','normal')
ylabel('Frequency (Hz)')
xlabel("Time (sec)")
xlim([t_min,t_max])
linkaxes(ax)

%%%% Moving pearson correlation
figure
set(gcf,'Units','Normalized','Position',[0 0.3 1 0.5])
set(gca,'Units','Normalized','Position',[0 0.1 1 0.8])

x = bstate(:,6);
[x_start,x_end,x_list] = sectionalize(x);
c = [0.6 0.2 0.2; % red: resting with the target bat
    0.3 0.7 0.9; % light blue: individual resting
    0 0 0; % flying
    0.1 0.6 0.2]; % green: resting with non-target bat
% c = [0.6 0.2 0.2; % ligh blue: resting with the target bat
%     0.3 0.7 0.9; % black: individual resting
%     0.0 0.0 0.0; % dark blue: flying
%     0.6 0.2 0.2; % green: resting with 1 neighbor
%     0.6 0.2 0.2]; % red: resting with 2>= neighbors
% Behavioral state
ax(1) = subplot(2,1,1);
hold on
for i = 1:length(x_list)
    xregion(t_lfp(x_start{i}),t_lfp(x_end{i}),'FaceColor',c(i,:),'FaceAlpha',0.5)
end
hold off
xlim([t_min t_max])
% LFP correlation
ax(2) = subplot(2,1,2);
hold on
plot(t_lfp,lfp_movcorr_1,'k')
plot(t_lfp,lfp_movcorr_2,'r')
hold off
xlim([t_min t_max])
ylabel('Amplitude')
xlabel("Time (sec)")
linkaxes(ax,'x')



%%%% proportion of variance
sz = 40;
figure;
subplot(1,2,1)
hold on
scatter(1,explained_1(1)/100,sz,'r','filled')
scatter(2,explained_1(2)/100,sz,'b','filled')
scatter(3:length(explained_1),explained_1(3:end)/100,sz,'k','filled')
hold off
xlabel('Principal component')
ylabel('Propotion of variance explained')
ylim([0 0.5])
title('Bat 1')

subplot(1,2,2)
hold on
scatter(1,explained_2(1)/100,sz,'r','filled')
scatter(2,explained_2(2)/100,sz,'b','filled')
scatter(3:length(explained_2),explained_2(3:end)/100,sz,'k','filled')
hold off
xlabel('Principal component')
ylabel('Propotion of variance explained')
ylim([0 0.5])
title('Bat 2')

%%%% Relative normalized power
figure; hold on
for i = 1:2
    plot(f_spctr_1,coeff_pca_1(:,i),'DisplayName',sprintf('PC%d',i));
end
hold off
legend()
ylabel('Relative normalized power')
xlabel('Frequency (Hz)')

%%%% PCA score
figure;
scatter(score_1(:,1),score_1(:,2))



%% Old scripts
% for i = 1:2
% 
%     figure('units','normalized','outerposition',[0 0.2 1 0.7]);
%     tiledlayout(3,1,'TileSpacing','tight')
%     ax(1) = nexttile; plot(t_corr/60,lfp_corr);
%     title(sprintf('LFP corr: Freq band = [%d %d], Date = %s, window = %dmin, lag = %dmin',freq_band_pc(i,1),freq_band_pc(i,2), batdate,wd/60,lag/60))
%     ylim([-0.3 0.3])
%     ylabel('Correlation')
%     ax(2) = nexttile;
%     % plot(t_1KHz_imp1(1e6:end)/60,LFP_imp1(1,1e6:end)); hold on
%     % plot(t_1KHz_imp2(1e6:end)/60,LFP_imp2(1,1e6:end))
%     plot(t_1KHz_imp1/60,LFP_imp1); hold on
%     plot(t_1KHz_imp2/60,LFP_imp2)
%     hold off
%     ylabel('LFP traces')
%     ylim([-4000 4000])
% 
%     ax(3) = nexttile;
%     plot(t_imp_resting/60,d_imp_resting)
%     xlabel('Time (min)')
%     ylabel('distance')
%     linkaxes(ax,'x')
% end