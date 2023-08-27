function RestPosClustering()

%% Load data
% Behavior
BehFile = dir('Ext_Behavior*\Extracted_Behavior*.mat');
load(fullfile(BehFile.folder,BehFile.name),'bflying','n_tags','T','r')

%% Clustering of resting position
% timepoint of landing and right before takeoff
tp_resting = cell(1,n_tags);
for bb = 1:n_tags
    tp_land = find(diff([bflying(1,bb);bflying(:,bb)]) == -1);
    tp_land = [1;tp_land];
    tp_toff = find(diff([bflying(:,bb);bflying(end,bb)]) == 1);
    tp_toff = [tp_toff;T];
    tp_resting{1,bb} = [tp_land tp_toff];
end

% mean position
for bb=1:n_tags
    for i=1:size(tp_resting{1,bb})
        tp_resting{1,bb}(i,3:5) = mean(r(tp_resting{1,bb}(i,1):tp_resting{1,bb}(i,2),:,bb));
    end
end

% clustering of mean positions concatenated across bats
r_mean_rest = [];
r_mean_rest_bid = [];
for bb = 1:n_tags
    r_mean_rest = [r_mean_rest; tp_resting{1,bb}(:,3:5)]; 
    r_mean_rest_bid = [r_mean_rest_bid; repmat(bb,[size(tp_resting{1,bb},1),1])];
end

% idx = kmeans(r_mean_rest,5);
idx = dbscan(r_mean_rest,0.25,5);

clu_unq = unique(idx);
figure
hold on
for cc = 1:length(clu_unq)
    plot3(r_mean_rest(idx==clu_unq(cc),1),r_mean_rest(idx==clu_unq(cc),2),r_mean_rest(idx==clu_unq(cc),3),'.')
    xlabel('x'); ylabel('y');zlabel('z')
    view([45 45])
    pause(0.5)
end
title('DBSCAN with epsilon 0.25 & minpts 5')
legend(num2str(clu_unq),'Location','northeastoutside')


for bb = 1:n_tags
    tp_resting{1,bb}(:,6) = idx(r_mean_rest_bid == bb);
end

mkdir(fullfile(BehFile.folder,'Analysis_Tatsumi'))
saveas(gcf,fullfile(BehFile.folder,'Analysis_Tatsumi','figure1.fig'))

%% Specify the reward cluster
prompt = {'Enter the reward cluster'};
clu_rwrd = inputdlg(prompt);
clu_rwrd = str2num(clu_rwrd{:});

% clu_rwrd = 2; % Cluster for reward position

%% Create output
r_clus_id = nan(size(bflying));
for bb = 1:n_tags
    for i = 1:size(tp_resting{1,bb},1)
        r_clus_id(tp_resting{1,bb}(i,1):tp_resting{1,bb}(i,2),bb) = tp_resting{1,bb}(i,6);
    end
end

%% Save
save(fullfile(BehFile.folder,'clu_pos.mat'),'r_clus_id','clu_rwrd')

%% Close file
close gcf