%% This is the log of scripts which is not used anymore

%% Manual clustering of landing sites
% Color setting for landing spots
bloc_clr = [0 0.4470 0.7410;...
    0.8500 0.3250 0.0980;...
    0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;...
    0.4660 0.6740 0.1880;...
    0.9098 0.1608 0.6235;...
    0.9686 0.1412 0.1137;
    0.0000 0.0000 0.0000];

%=== FIGURE: check the landing sites on data
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0.2 0.25 0.7 0.35]);
tiledlayout(1, 4, 'TileSpacing', 'tight');
% Top view
nexttile;hold on; for i=1:n_tags
    scatter3(r(logical(1-bflying(:,i)),1,i),r(logical(1-bflying(:,i)),2,i),r(logical(1-bflying(:,i)),3,i),'.','MarkerEdgeColor',bat_clr(i,:));
    xlabel("x"); ylabel("y"); view([90 90]);
end; hold off
% Side view 1
nexttile; hold on; for i=1:n_tags
    scatter3(r(logical(1-bflying(:,i)),1,i),r(logical(1-bflying(:,i)),2,i),r(logical(1-bflying(:,i)),3,i),'.','MarkerEdgeColor',bat_clr(i,:));
    xlabel("x"); ylabel("y"); zlabel("z"); view([0 0]);
end; hold off
% Side view 2
nexttile; hold on; for i=1:n_tags
    scatter3(r(logical(1-bflying(:,i)),1,i),r(logical(1-bflying(:,i)),2,i),r(logical(1-bflying(:,i)),3,i),'.','MarkerEdgeColor',bat_clr(i,:));
    xlabel("x"); ylabel("y"); zlabel("z"); view([90 0]);
end; hold off
% Corner view
nexttile; hold on; for i=1:n_tags
    scatter3(r(logical(1-bflying(:,i)),1,i),r(logical(1-bflying(:,i)),2,i),r(logical(1-bflying(:,i)),3,i),'.','MarkerEdgeColor',bat_clr(i,:));
    xlabel("x"); ylabel("y"); zlabel("z"); view([-45 45]);
end; hold off

% Ad-hoc allocation of edges
bloc_edge = {[-3 -1.5], [1.5 3], [1.25 2.5];... % site 1
    [1.8 3], [1.5 3], [1.25 2.5];... % site 2
    [1.8 3], [0.5 1.5], [1.25 2.5];... % site 3
    [1.8 3], [1 2], [0.25 1.25];... % site 4
    [1.8 3], [-1.5 -0.5], [1.25 2];... % site 5
    [1.8 3], [-2.5 -1], [0.25 1];... % site 6
    [1.8 3], [-3 -1.5], [1.5 2.5]}; % site 7

nloc = size(bloc_edge,1)+1; % number of landing sites

% Color setting for landing sites
bloc_clr = [0 0.4470 0.7410;...
    0.8500 0.3250 0.0980;...
    0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;...
    0.4660 0.6740 0.1880;...
    0.9098 0.1608 0.6235;...
    0.9686 0.1412 0.1137;
    0.0000 0.0000 0.0000];

% Ad-hoc segmentation: 0=flying, 1~nloc-1=landing sites, nloc=others
bloc = zeros(T,n_tags);
for i=1:n_tags
    bloc(logical(1-bflying(:,i)),i) = 8; % cluster 8 (unusual positions)
    for j=1:size(bloc_edge,1)
        bloc(logical(1-bflying(:,i))& ...
            r(:,1,i)>bloc_edge{j,1}(1) & r(:,1,i)<bloc_edge{j,1}(2) & ... % x
            r(:,2,i)>bloc_edge{j,2}(1) & r(:,2,i)<bloc_edge{j,2}(2) & ... % y
            r(:,3,i)>bloc_edge{j,3}(1) & r(:,3,i)<bloc_edge{j,3}(2), i) = j;  % z
        bloc_cnt(j,i) = sum(bloc(:,i)==j);
    end
end

% Probability of landing at each site and flying
bloc_cnt = zeros(nloc+1,n_tags); bloc_cnt(1,:)=sum(bflying);
for i=1:nloc
    bloc_cnt(i+1,:) = sum(bloc(:,:)==i);
end
bloc_prb = bloc_cnt / T; % probability

%=== FIGURE: check the segmentation
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0.2 0.25 0.7 0.35]);
tiledlayout(1, 4, 'TileSpacing', 'tight');
% Top view
nexttile; hold on
for i=1:n_tags
    for j=1:nloc+1
        if sum(bloc(:,i)==j)>0
            scatter3(r(bloc(:,i)==j,1,i),r(bloc(:,i)==j,2,i),r(bloc(:,i)==j,3,i),'.','MarkerEdgeColor',bloc_clr(j,:));
            xlabel("x"); ylabel("y"); zlabel("z"); view([90 90]);
        end
    end
end; hold off
% Side view 1
nexttile; hold on
for i=1:n_tags
    for j=1:nloc+1
        if sum(bloc(:,i)==j)>0
            scatter3(r(bloc(:,i)==j,1,i),r(bloc(:,i)==j,2,i),r(bloc(:,i)==j,3,i),'.','MarkerEdgeColor',bloc_clr(j,:));
            xlabel("x"); ylabel("y"); zlabel("z"); view([0 0]);
        end
    end
end; hold off
% Side view 2
nexttile; hold on
for i=1:n_tags
    for j=1:nloc+1
        if sum(bloc(:,i)==j)>0
            scatter3(r(bloc(:,i)==j,1,i),r(bloc(:,i)==j,2,i),r(bloc(:,i)==j,3,i),'.','MarkerEdgeColor',bloc_clr(j,:));
            xlabel("x"); ylabel("y"); zlabel("z"); view([90 0]);
        end
    end
end; hold off
% Corner view
nexttile; hold on
for i=1:n_tags
    for j=1:nloc+1
        if sum(bloc(:,i)==j)>0
            scatter3(r(bloc(:,i)==j,1,i),r(bloc(:,i)==j,2,i),r(bloc(:,i)==j,3,i),'.','MarkerEdgeColor',bloc_clr(j,:));
            xlabel("x"); ylabel("y"); zlabel("z"); view([-45 45]);
        end
    end
end; hold off

%=== FIGURE: visualizing the site preference
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0.3 0.25 0.5 0.35]);
bar(bloc_prb(2:nloc,:)'); ylabel("Probability"); xticklabels(bat_nms); legend({'1','2','3','4','5','6','7'})

%== FIGURE: hierarchical clustering of landing spots
figure();dendrogram(linkage(pdist(bloc','hamming')));xlabel("bat ID")
figure();dendrogram(linkage(pdist(bloc_cnt')));xlabel("bat ID")

%% k-means clustsering
% Specification of initial positions of centroids
n_spots = 7;
cent_ini = [-2.2 2.2 2;...
    2.8 2.5 2.1;...
    2.2 -2.5 2.1;...
    2.8 -1 1.7;...
    2.8 1 1.7;...
    2.8 -1.5 0.8;...
    2.8 1.4 0.8];

% concatenate all bats
xy = reshape(permute(r,[1 3 2]),[T*n_tags 3]); 

%=== clustering
tic; [cid_a,C] = kmeans(xy(:,:),n_spots,'Start',cent_ini); % 1~: cluster; no detection of outliers
cid_a(logical(reshape(bflying,[T*n_tags 1]))) = -1; % Exclude the positions during flying
cid_a(xy(:,1)<1.7 & xy(:,2)<1.51) = -1; % Exclude the start positions
idx = reshape(cid_a, [T n_tags]); % cluster index [time, bats]

%=== FIGURE: 3D visualization of landing spot clusters
figure(); hold on
for j=1:size(C,1)
    plot3(xy(cid_a==j,1),xy(cid_a==j,2),xy(cid_a==j,3),'.')
    xlabel("x"); ylabel("y"); zlabel("z"); view([135 45]);
end
plot3(xy(idx==-1,1),xy(idx==-1,2),xy(idx==-1,3),'k.');
hold off;
toc; fprintf('%d cluster was identified! \n', length(unique(cid_a))-1);

save('kmeans_idx.mat','idx','C')

%% DBSCAN for 1 bat
tic
i=1; n=300000;
[idx C] = dbscan(r(1:n,:,i),0.13,100);
idx(logical(bflying(1:n,i))) = -1; % Exclude the positions during flying
idx(r(1:n,1,i)<1.7 & r(1:n,2,i)<1.51) = -1; % Exclude the start positions
figure(); hold on
for j=0:max(idx)
%     scatter3(r_land(idx==j,1),r_land(idx==j,2),r_land(idx==j,3));
    plot3(r(idx==j,1),r(idx==j,2),r(idx==j,3),'.')
    xlabel("x"); ylabel("y"); zlabel("z"); view([135 45]);
end
plot3(r(idx==-1,1),r(idx==-1,2),r(idx==-1,3),'k.')
hold off;
toc; fprintf('%d cluster was identified! \n', length(unique(idx))-1);

%% DBSCAN for all bats
% tic
% % concatenate all bats
% xy = reshape(permute(r,[1 3 2]),[T*n_tags 3]); 
% %=== Clustering
% [idx_a C] = dbscan(xy(:,:),0.13,100);
% idx_a(logical(reshape(bflying,[T*n_tags 1]))) = -1; % Exclude the positions during flying
% idx_a(xy(:,1)<1.7 & xy(:,2)<1.51) = -1; % Exclude the start positions
% idx = reshape(idx_a, [T n_tags]); % cluster index [time, bats]
% 
% %=== FIGURE: 3D visualization of landing spot clusters
% figure(); hold on
% for j=1:size(C,1)
%     plot3(xy(idx_a==j,1),xy(idx_a==j,2),xy(idx_a==j,3),'.')
%     xlabel("x"); ylabel("y"); zlabel("z"); view([135 45]);
% end
% % plot3(xy(idx==-1,1),xy(idx==-1,2),xy(idx==-1,3),'k.');
% hold off; toc;
% fprintf('%d cluster was identified! \n', length(unique(idx_a))-1);
% 
% save('dbscan_idx.mat', 'idx', 'C')

%% Fast DBSCAN using kdtrees
% tic
% i=1;
% idx = fdbscan(r(1:300000,:,i),10,0.1);
% idx(logical(bflying(:,i)),:,i) = 0;
% %==FIGURE
% figure(); hold on
% for j=1:max(idx)
% %     scatter3(r_land(idx==j,1),r_land(idx==j,2),r_land(idx==j,3));
%     plot3(r(idx==j,1),r(idx==j,2),r(idx==j,3),'.')
%     xlabel("x"); ylabel("y"); zlabel("z"); view([-45 45]);
% end
% plot3(r(idx==0,1),r(idx==0,2),r(idx==0,3),'k.')
% hold off; toc; fprintf('%d cluster was identified!', length(unique(idx))-1)

%% DBSCAN with knnsearch


eval_cols = ["whole data size", "epsilon", "minpts",...
    "# of clusters", "train outliers (%)", "all outliers (%)", "processing time"];
eval_clustering = zeros(4,7);
eval_clustering(:,1) = sum(1-bflying(:));
eval_clustering(:,2) = [0.10:0.15];
eval_clustering(:,3) = eval_clustering(:,2)./eval_clustering(:,1);

for i=1:size(eval_clustering,1)
    fprintf('Clustering for %d points is processing... \n', eval_clustering(i,2))
    tic
    % Concatenate all sessions
    r_all = reshape(permute(r,[1 3 2]),[T*n_tags 3]); % concatenated positions for all bats
    idx_all = reshape(repmat(1:n_tags,T,1),[T*n_tags 1]); % concatenated index for all bats
    bflying_all = reshape(bflying, [T*n_tags 1]); % concatenated flying states for all bats

    % Downsample concatenated data
    idx_train = sort(randsample(find(1-bflying_all),eval_clustering(i,2))); % index of positions for clustering
    idx_test = find(1-bflying_all); idx_test = idx_test(~ismember(idx_test,idx_train)); % index of positions for knn-search

    % DBSCAN for downsampled data
    [cid_train C] = dbscan(r_all(idx_train,:),0.13,100); % dbscan clusteirng 
    fprintf('%d cluster was identified! \n', length(unique(cid_train(cid_train~=-1))));
    % knnsearch for other session
    cid_test = cid_train(knnsearch(r_all(idx_train,:), r_all(idx_test,:)));

    % concatenate train and test data
    [idx_brest, sortidx] = sort([idx_train;idx_test]); cid_all = [cid_train; cid_test]; cid_all = cid_all(sortidx);
    
    % evaluation
    eval_clustering(i,4) = length(unique(cid_train(cid_train~=-1)));
    eval_clustering(i,5) = sum(cid_train == -1)/eval_clustering(i,2);
    eval_clustering(i,6) = sum(cid_all == -1)/eval_clustering(i,1);
    eval_clustering(i,7) = toc;

    %=== FIGURE: 3D visualization of landing spot clusters
    figure(); hold on
    for j=1:length(unique(cid_all))-1
        plot3(r_all(idx_brest(cid_all==j),1),r_all(idx_brest(cid_all==j),2),r_all(idx_brest(cid_all==j),3),'.')
    end
    plot3(r_all(idx_brest(cid_all==-1),1),r_all(idx_brest(cid_all==-1),2),r_all(idx_brest(cid_all==-1),3),'k.');
    xlabel("x"); ylabel("y"); zlabel("z"); view([135 45]);
    hold off;

    %=== Save
    saveas(gcf,strcat('Figs_Tatsumi/knn_dbscan_',num2str(eval_clustering(i,2))),'png')
end

%% Spectral clustering
% tic
% i=1;
% idx = spectralcluster(r(1:30000,:,i),10 );
% % idx(logical(bflying(:,i)),:,i) = 0;
% %==FIGURE
% figure(); hold on
% for j=0:max(idx)
% %     scatter3(r_land(idx==j,1),r_land(idx==j,2),r_land(idx==j,3));
%     plot3(r(idx==j,1),r(idx==j,2),r(idx==j,3),'.')
%     xlabel("x"); ylabel("y"); zlabel("z"); view([-45 45]);
% end
% plot3(r(idx==0,1),r(idx==0,2),r(idx==0,3),'k.')
% hold off; toc; fprintf('%d cluster was identified!', length(unique(idx))-1)

%% knn search for each session
r_all = reshape(permute(r,[1 3 2]),[T*n_tags 3]); % concatenated positions for all bats
bflying_all = reshape(bflying, [T*n_tags 1]); % concatenated flying states for all bats

r_brest = r_all(logical(1-bflying_all),:); % positions at rest

% knn search
cid = cid_temp(knnsearch(r_temp(:,:), r_brest(:,:)));
    
%=== FIGURE: 3D visualization of landing spot clusters
figure(); hold on
xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:)); xlabel("x"); ylabel("y"); zlabel("z"); view([135 45]);
for j=1:length(unique(cid))-1
    plot3(r_brest(cid==j,1),r_brest(cid==j,2),r_brest(cid==j,3),'.')
end
plot3(r_brest(cid==-1,1),r_brest(cid==-1,2),r_brest(cid==-1,3),'k.');

hold off;


%% Labeling of behavioral states
% cluster 1-3: home, cluster 4-7: feeders
% idx_home = idx <= 3 & idx ~= -1; % 1=home
% idx_feed = idx >= 4; % 1=feeders
% idx_soc = false(size(idx)); % 1=social (more than one bat for one place), 0=unsocial
% for i=1:T
%     for j=1:n_tags
%        if sum(idx(i,:)==idx(i,j)) >= 2 && idx(i,j)~=-1 % avoid the case bats rests at abnormal places or flying
%            idx_soc(i,j) = true;
%        end
%     end
% end
% 
% %=== labeling
% % 1=home-social, 2=home-unsocial, 3=feeder-social, 4=feeder-unsocial, 5=flying OR abnormal resting spots
% idx_state = [1, 2, 3, 4, 5];
% num_state = length(idx_state);
% states = {'others';'home-social';'home-unsocial';'feeder-social';'feeder-unsocial'};
% bstate = zeros(size(idx));
% bstate(idx~=-1&idx_home&idx_soc) = idx_state(1);
% bstate(idx~=-1&idx_home&~idx_soc) = idx_state(2);
% bstate(idx~=-1&idx_feed&idx_soc) = idx_state(3);
% bstate(idx~=-1&idx_feed&~idx_soc) = idx_state(4);
% bstate(idx==-1) = idx_state(5);


%% Count transitions of takeoff -> landing state
statenames = ["social-home" "isolated-home" "social-feeder" "isolated-feeder"];
blanding = false(T,n_tags); % true=landing
btakeoff = false(T,n_tags); % true=takeoff

P = zeros(4,4,n_tags); % transition matrix
bloc_cnt = zeros(length(unique(idx))-1,n_tags);
for i=1:n_tags
    blanding(:,i) = diff([-1;idx(:,i)])~=0 & idx(:,i)~=-1; % logical vector for the i-th bat's landing timing
    bspots = idx(blanding(:,i),:); % spots transition vector for the i-th bat
    for j=1:size(bloc_cnt,1)
        bloc_cnt(j,i) = sum(bspots(:,i)==j);
    end
    bstate = zeros(size(bspots,1),1); % state transition vector for the i-th bat
    for j=1:length(bstate)
        bspot_now = bspots(j,i); % current spot
        % label current state
        if bspot_now<=3 && bspot_now~=-1 && sum(bspots(j,:)==bspot_now) >= 2
            bstate(j) = 1; % social-home: home (1-3) and social (>= 2 animals)
        elseif bspot_now<=3 && bspot_now~=-1 && sum(bspots(j,:)==bspot_now)==1
            bstate(j) = 2; % unsocial-home: home (1-3) and unsocial (< 2 animlas)
        elseif bspot_now>=4 && sum(bspots(j,:)==bspot_now)>=2
            bstate(j) = 3; % social-feeder: feeder (4-7) and social (>= 2 animals)
        elseif bspot_now>=4 && sum(bspots(j,:)==bspot_now)==1
            bstate(j) = 4; % unsocial-feeder: feeder (4-7) and unsocial (< 2 animals)
        end
    end

    % Count state transitions
    for j=1:length(bstate)-1
        P(bstate(j),bstate(j+1),i) = P(bstate(j),bstate(j+1),i) + 1;
    end
end

%=== FIGURE: visualizing the site preference
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0.3 0.25 0.5 0.35]);
bar(bloc_cnt(:,:)'); ylabel("Probability"); xticklabels(bat_nms); legend({'Home 1','Home 2','Home 3','Feeder 4','Feeder 5','Feeder 6','Feeder 7'},'Location','northeastoutside')

%=== FIGURE: visualization of state transition
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0.05 0.1 0.9 0.8]);
tiledlayout(2, ceil(n_tags/2), 'TileSpacing', 'tight');
for i=1:n_tags
    mc = dtmc(P(:,:,i),'StateNames',statenames);
    nexttile(); graphplot(mc,'ColorEdges',true)
    title(bat_nms(i,:),'FontSize',12,'Color',bat_clr(i,:))
end

%=== FIGURE: visualization of state transition
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0.05 0.2 0.9 0.3]);
tiledlayout(1, n_tags);
for i=1:n_tags
    nexttile(); imagesc(P(:,:,i)); colormap(jet); title(bat_nms(i,:),'FontSize',12,'Color',bat_clr(i,:))
    yticks([1:4])
    if i>1
        yticks([]);
    end
end

%% Count State Transitions for takeoff
% statenames = ["social-home" "isolated-home" "social-feeder" "isolated-feeder"];
% btakeoff = false(T,n_tags);
% P = zeros(4,4,n_tags);
% for i=1:n_tags
%     btakeoff(:,i) = diff([idx(:,i);0])~=0 & idx(:,i)~=-1; % logical vector for the i-th bat's landing timing
%     btakeoff(end,i) = false;
%     bspots = idx(btakeoff(:,i),:); % spots transition vector for the i-th bat
%     bstate = zeros(size(bspots,1),1); % state transition vector for the i-th bat
%     for j=1:length(bstate)
%         bspot_now = bspots(j,i); % current spot
%         % label current state
%         if bspot_now<=3 && bspot_now~=-1 && sum(bspots(j,:)==bspot_now) >= 2
%             bstate(j) = 1; % social-home: home (1-3) and social (>= 2 animals)
%         elseif bspot_now<=3 && bspot_now~=-1 && sum(bspots(j,:)==bspot_now)==1
%             bstate(j) = 2; % unsocial-home: home (1-3) and unsocial (< 2 animlas)
%         elseif bspot_now>=4 && sum(bspots(j,:)==bspot_now)>=2
%             bstate(j) = 3; % social-feeder: feeder (4-7) and social (>= 2 animals)
%         elseif bspot_now>=4 && sum(bspots(j,:)==bspot_now)==1
%             bstate(j) = 4; % unsocial-feeder: feeder (4-7) and unsocial (< 2 animals)
%         end
%     end
% 
%     % Count state transitions
%     for j=1:length(bstate)-1
%         P(bstate(j),bstate(j+1),i) = P(bstate(j),bstate(j+1),i) + 1;
%     end
% end
% 
% %=== FIGURE: visualization of state transition
% figure(); set(gcf, 'units', 'normalized', 'outerposition', [0.05 0.1 0.9 0.8]);
% tiledlayout(2, ceil(n_tags/2), 'TileSpacing', 'tight');
% for i=1:n_tags
%     mc = dtmc(P(:,:,i),'StateNames',statenames);
%     nexttile(); graphplot(mc,'ColorEdges',true)
% end

%%
% Extract the timings of landings and takeoffs
% [btakeoff,blanding] = findchanges(bposnew,[-1 0],1);
% 
% P = zeros(4,4,n_tags); % transition matrix
% 
% for i=1:n_tags    
%     Identify the spots when bats takeoff and land
%     bspot_land = bposnew(blanding(:,i),:); 
%     bspot_land = bspot_land(2:end,:); % remove the landing to initial spot
%     bspot_toff = bposnew(btakeoff(:,i),:);
%     if ismember(bposnew(end,i),[-1 0]) % if the i-th bat flied when task finished
%         bspot_toff = bspot_toff(1:end-1); % remove the final takeoff
%     end
%     
%     bspots = zeros(size(bspot_land,1),2,n_tags); % column 1 is takeoff spot, column 2 is landing spot for D2
%     bspots(:,1,:) = bspot_toff; 
%     bspots(:,2,:) = bspot_land;
% 
%     Label the takeoff and landing spots
%     bstate = zeros(size(bspots,1),2); % state transition vector for the i-th bat
%     for j=1:length(bstate)
%         for k=1:2 % takeoff and landing
%             bspot_now = bspots(j,k,i); % current spot (takeoff/landing)
%             label current state
%             if ismember(bspot_now,11:19) && sum(bspots(j,k,:)==bspot_now)>=2
%                 bstate(j,k) = 1; % social-home: home (1-3) and social (>= 2 animals)
%             elseif ismember(bspot_now,11:19) && sum(bspots(j,k,:)==bspot_now)==1
%                 bstate(j,k) = 2; % unsocial-home: home (1-3) and unsocial (< 2 animlas)
%             elseif ismember(bspot_now,21:29) && sum(bspots(j,k,:)==bspot_now)>=2
%                 bstate(j,k) = 3; % social-feeder: feeder (4-7) and social (>= 2 animals)
%             elseif ismember(bspot_now,21:29) && sum(bspots(j,k,:)==bspot_now)==1
%                 bstate(j,k) = 4; % unsocial-feeder: feeder (4-7) and unsocial (< 2 animals)
%             end
%         end
%     end
% 
%     % Count state transitions
%     for j=1:size(bstate,1)
%         P(bstate(j,1),bstate(j,2),i) = P(bstate(j,1),bstate(j,2),i) + 1;
%     end
% end