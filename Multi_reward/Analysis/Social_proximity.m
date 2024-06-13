recID = [3,6,9,12,15,19];
p = nan(10,10,length(recID));
proximity_corr_1 = zeros(length(recID),1);
proximity_corr_2 = zeros(length(recID)-1,1);
for i = 1:length(recID)
    rootdir = fullfile(ParDir,allRecs.rootdir{recID(i)}); % root directory of the recording
    load(fullfile(rootdir,'analysis','social_proximity.mat'),'social_proximity','p_social_proximity','proximity_index')
    p(:,:,i) = proximity_index;
end

%% Corr day vs mean
% p_mean = mean(p,3);
for i = 1:length(recID)
    x = reshape(p(:,:,i),[1,100]); x= x(~isnan(x));
    p_mean = mean(p(:,:,[1:i-1,i+1:end]),3);
    y = reshape(p_mean,[1,100]); y = y(~isnan(y));
    xy_corr = corrcoef(x,y);
    proximity_corr_1(i) = xy_corr(1,2);
end

%%% FIGURE: proximity index correlation (each day vs mean)
fig = figure;
plot(1:length(recID),proximity_corr_1,'k','LineWidth',2)
set(gca,'FontSize',16)
xlabel('Day')
ylabel('Proximity Index Corr')
xlim([0.7 length(recID)+0.3])
ylim([-0.5 0.5])

saveas(fig,'Z:\users\Tatsumi\data\MultiReward\Analysis\figure\Social_proximity\social_proximity_index_corr_vsMean.jpg')

%% Corr day i vs day i+1
for i = 2:length(recID)
    x = reshape(p(:,:,i-1),[1,100]); x= x(~isnan(x));
    y = reshape(p(:,:,i),[1,100]); y= y(~isnan(y));
    xy_corr = corrcoef(x,y);
    proximity_corr_2(i-1) = xy_corr(1,2);
end

%%% FIGURE: proximity index correlation (each day vs mean)
fig = figure;
plot(2:length(recID),proximity_corr_2,'k','LineWidth',2)
set(gca,'FontSize',16)
xlabel('Day')
ylabel('Proximity Index Corr')
xlim([1.7 length(recID)+0.3])
ylim([-0.5 0.5])

saveas(fig,'Z:\users\Tatsumi\data\MultiReward\Analysis\figure\Social_proximity\social_proximity_index_corr_vsPreviousDay.jpg')

%%
%%% FIGURE
fig = figure();
fontsize(fig, 20, "points")
set(fig, 'units','normalized','Position',[0.2 0.2 0.55 0.55]);
tiledlayout(2,3,"TileSpacing",'compact');
for i = 1:length(recID)
    nexttile
    C = redblue(256);
    h = imagesc(p(:,:,i)/max(p(:,:,i),[],'all'));
    set(gca,'FontSize',16)
    set(h, 'AlphaData', ~isnan(proximity_index))
    clim([-1 1])
    colorbar
    colormap(C)
    % xlabel('Bat')
    % ylabel('Bat')
    title(sprintf('Day %d',i))
    xticks(1:10)
    yticks(1:10)
end

saveas(fig,'Z:\users\Tatsumi\data\MultiReward\Analysis\figure\Social_proximity\Social_proximity_progression.jpg')