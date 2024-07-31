% Population coding

%% Load data
recID = 12;
dataTypes = 'barsm'; % data types to be loaded
[BhvData,AnalogData,RewardData,spikeData,mtData] = LoadExperimentData(recID,dataTypes);
[Flights, TrjCluster] = extractFlights(BhvData,RewardData,mtData,'include','all');

%% Calculate the correlation between units
t_bhv_max = BhvData.t(end);
t_ephys_max = max(spikeData.allGlobalSpikeTimes_usec)/ 1e6;

t_min = 0;
t_max = min([t_bhv_max,t_ephys_max]);

bin = 0.1;
wd = 0.5;
t_bin = t_min:bin:t_max; % time bin for neural data

neurons = find(strcmp(spikeData.good_units.region,'MEC'));
d_corr = NaN(length(neurons),length(neurons));
spcnt = NaN(length(neurons),length(t_bin)-1);
frs = NaN(size(spcnt));
for i = 1:length(neurons)
    spcnt(i,:) = histcounts(spikeData.good_units.spikeTimes_usec{neurons(i)}/1e6,t_bin);
    frs(i,:) = movsum(spcnt(i,:),wd/bin,2)/0.5;
end
for i = 1:length(neurons)
    for j = 1:length(neurons)
        if i~=j
            corr_now = corrcoef(frs(i,:),frs(j,:));
            d_corr(i,j) = corr_now(1,2);
        else
            d_corr(i,j) = 0;
        end
    end
end

%% k-means clustering
rng(123); % For reproducibility
k = 2;
idx_kmean = kmeans(d_corr,k);
[idx_sorted,ia] = sort(idx_kmean);
d_corr_new = d_corr(ia,:); d_corr_new = d_corr_new(:,ia);
frs_sorted = frs(ia,:);

figure
set(gcf,'Units','normalized','Position',[0.3 0.3 0.45 0.3])
subplot(1,2,1)
imagesc(d_corr)
title('Raw')
xlabel('Neurons')
ylabel('Neurons')

subplot(1,2,2)
imagesc(d_corr_new)
title('Clustered')
xlabel('Neurons')
ylabel('Neurons')


% Sorted firing rate traces
figure
set(gcf,'Units','normalized','Position',[0 0.3 1 0.3])
imagesc(frs)
ylabel('Zscored firing rate (Hz)')
xlabel('Time')
set(gca,'FontSize',16)

%% Average firing activity for each ensemble
fr_mean = NaN(k,length(t_bin)-1);
for i = 1:k
    fr_mean(i,:) = normalize(mean(frs(idx_kmean==i,:),1),'zscore');
end

figure;
hold on

plot(t_bin(1:end-1)+bin/2,fr_mean(1,:),'b')
plot(t_bin(1:end-1)+bin/2,fr_mean(2,:),'r')
% plot(t_bin(1:end-1)+bin/2,fr_mean(3,:),'k')
hold off


%% Define behavioral state
n_bats = BhvData.n_tags;
T = BhvData.T;
r = BhvData.r;
bflying = BhvData.bflying;
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
    % Resting with other bats
    bstate(blanding(:,bb)==1 & n_neighbors(:,bb)>0,bb) = 1; % social resting
    % Resting alone
    bstate(blanding(:,bb)==1 & n_neighbors(:,bb)==0,bb) = 2; % # of neighbors
end

%% Align behavioral state to neural recording
t = BhvData.t;
n_bats = BhvData.n_tags;
bstate_new = NaN(length(t_bin)-1,n_bats);

for bb = 1:n_bats
    bstate_new(:,bb) = interp1(t,bstate(:,bb),t_bin(1:end-1)+bin/2,'nearest');
end

bstate = bstate_new;



%% FIGURE
t_fig = t_bin(1:end-1)+bin/2;
figure
set(gcf,'Units','Normalized','Position',[0 0.3 1 0.5])
set(gca,'Units','Normalized','Position',[0 0.1 1 0.8])
% Behavioral state
x = bstate(:,6);
[x_start,x_end,x_list] = sectionalize(x);
c = [0 0 0; % black: flying
    0.6 0.2 0.2; % red: resting with the target bat
    0.3 0.7 0.9]; % light blue: individual resting

ax(1) = subplot(2,1,1);
hold on
for i = 1:length(x_list)
% for i = 4
    xregion(t_fig(x_start{i}),t_fig(x_end{i}),'FaceColor',c(i,:),'FaceAlpha',0.5)
end
hold off
xlim([t_min t_max])
yticks([])
% LFP spectrogram
ax(2) = subplot(2,1,2);
imagesc(frs)
ylabel('Zscored firing rate (Hz)')
xlabel('Time')
set(gca,'FontSize',16)
xlim([t_min t_max])

linkaxes(ax)