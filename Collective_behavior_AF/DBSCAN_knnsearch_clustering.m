function DBSCAN_knnsearch_clustering(opt,overwrite)
% This script is for creating a template to classify bats' positions.
% Data of positions for all bats in all sessions are concatenated and
% downsampled. The downsampled data is then classified into several
% clusters by DBSCAN. The cluster identification of these positions can be
% used for determining the cluster identification of unlabeled data by
% k-nearest neighbor method
%
% Inputs
%     opt: parameters for DBSCAN
%     overwrite: if overwrite is true, then clustering_template.mat would
%     be updated
% Output
%     cid_temp: cluster of positions
%     C_temp: centroid positions
%     r_temp: positions
%     epsilon: an epsilon estimated by knee method
%
% Note 04/05/2023 by TY
% bout.mat 

smpts = opt.smpts; % sampling size; % parameter of dbscan
minpts = opt.minpts; % parameter of dbscan
th = opt.th; % threshold for exclusion of small clusters

sessions = dir('**/Extracted_Behavior*.mat'); % extract file paths

%% DBSCAN
if exist(fullfile(cd,'clustering_template.mat'),'file') && ~overwrite
    load(fullfile(cd,'clustering_template.mat'))
else
    % Concatenate all sessions
    for i=1:length(sessions)
        load(fullfile(sessions(i).folder, sessions(i).name),'r','T','n_tags','bflying','r_lim') % load each session
        r_session = reshape(permute(r,[1 3 2]),[T*n_tags 3]); % align all bats' positions
        bflying_session = reshape(bflying, [T*n_tags 1]); % align all bats' bflying

        if i==1
            r_all = r_session;
            bflying_all = bflying_session;
        else
            r_all = [r_all; r_session]; % concatenate all sessions
            bflying_all = [bflying_all; bflying_session]; % concatenate all sessions
        end
    end

    % downsampling
    idx = sort(randsample(find(1-bflying_all),smpts)); % index of positions for clustering
    r_temp = r_all(idx,:); % downsampled positions

    % Find the optimal epsilon
    kD = pdist2(r_temp,r_temp,'euclidean','Smallest',minpts);
    Y = sort(kD(end,:));
    X_knee = knee_pt(Y); % knee-pt
    % X_knee = findchangepts(Y,'Statistic','linear'); % Find the optimal epsilon by findchangepts 
    epsilon = Y(X_knee); 

    %=== FIGURE* k-dist plot
    figure(); plot(Y); hold on; plot(X_knee, epsilon,'ro'); hold off
    title('k-distance graph'); xlabel('Points sorted with nearest distances'); ylabel('Nearest distances')
    saveas(gcf,fullfile(cd,'kdist.png'))

    % DBSCAN clustering
    fprintf('Processing DBSCAN with epsilon %.2f, minpts %d smpts %d \n', epsilon, minpts, smpts);
    [cid_temp, C_temp] = dbscan(r_temp,epsilon,minpts); % dbscan clusteirng 
    fprintf('%d clusters were detected! \n', length(unique(cid_temp(cid_temp~=-1))));

    % Exclude small clusters
    fprintf('Excluding small clusters... \n')
    n_rm = 0; % count of removed clusters
    new_id = 1; % new id for remained clusters
    for j=1:length(unique(cid_temp))-1
       if sum(cid_temp == j) < length(cid_temp)*th
           cid_temp(cid_temp==j) = -1; % remove small clusters 
           n_rm = n_rm + 1;
       else
           cid_temp(cid_temp == j) = new_id; % reallocate cluster numbers
           new_id = new_id + 1;
       end
    end
    fprintf('%d clusters were removed. \n', n_rm)
    fprintf('%d clusters was finally identified. \n', length(unique(cid_temp))-1)


    %=== FIGURE: 3D visualization of landing spot clusters
    figure(); hold on
    box on; ax = gca; ax.BoxStyle = 'full';
    xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:)); xlabel("x"); ylabel("y"); zlabel("z"); view([135 50]);
    title(sprintf('DBSCAN: eps=%.2f, minpts=%d, smpts=%dk',epsilon,minpts,smpts/1000))
    for k=1:length(unique(cid_temp))-1
    %     plot3(r_temp(cid_temp==k,1),r_temp(cid_temp==k,2),r_temp(cid_temp==k,3),'.')
        scatter3(r_temp(cid_temp==k,1),r_temp(cid_temp==k,2),r_temp(cid_temp==k,3),'.','MarkerEdgeAlpha',0.65)
    %     scatter3(r_temp(cid_temp==k,1),r_temp(cid_temp==k,2),r_temp(cid_temp==k,3),5,'filled','MarkerFaceAlpha',0.2)
        saveas(gcf,fullfile(cd,['rest_cluster_temp_' num2str(k) '.png']))
    end
    % plot3(r_temp(cid_temp==-1,1),r_temp(cid_temp==-1,2),r_temp(cid_temp==-1,3),'k.');
    % scatter3(r_temp(cid_temp==-1,1),r_temp(cid_temp==-1,2),r_temp(cid_temp==-1,3),'k.','MarkerEdgeAlpha',0.2);
    % scatter3(r_temp(cid_temp==-1,1),r_temp(cid_temp==-1,2),r_temp(cid_temp==-1,3),'.','MarkerEdgeColor','k','MarkerEdgeAlpha',0.03);
    scatter3(r_temp(cid_temp==-1,1),r_temp(cid_temp==-1,2),r_temp(cid_temp==-1,3),5,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.03);

    hold off;
    saveas(gcf,fullfile(cd,'rest_cluster_temp.png'))

    %=== Save
    save(fullfile(cd,'clustering_template.mat'),'r_temp','cid_temp','C_temp','epsilon','smpts','minpts','th') % template

end

%% knn-search
% position clustering by knn-search
for i=1:length(sessions)
    fprintf('Current session is %s \n', sessions(i).folder)
    load(fullfile(sessions(i).folder, sessions(i).name),'r','T','n_tags','bflying','r_lim','bat_nms','bat_clr') % load each session

    if exist(fullfile(sessions(i).folder,'bpos.mat'),'file') && ~overwrite
        load(fullfile(sessions(i).folder,'bpos.mat'))
    else
    
        % knn search for each session
        r_all = reshape(permute(r,[1 3 2]),[T*n_tags 3]); % concatenated positions for all bats
        bflying_all = reshape(bflying, [T*n_tags 1]); % concatenated flying states for all bats
        r_brest = r_all(logical(1-bflying_all),:); % positions at rest

        cid = cid_temp(knnsearch(r_temp(:,:), r_brest(:,:))); % knn search

        % Merge with bflying
        bpos = 1- bflying_all; % bpos: 0 = flying, -1,1-k = cluster ID of resting spots
        bpos(bpos~=0) = cid;
        bpos = reshape(bpos, [T n_tags]);

        %=== FIGURE: 3D visualization of landing spot clusters
        figure(); hold on
        box on; ax = gca; ax.BoxStyle = 'full';
        xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:)); xlabel("x"); ylabel("y"); zlabel("z"); view([135 50]);
        for j=1:length(unique(cid))-1
    %         plot3(r_brest(cid==j,1),r_brest(cid==j,2),r_brest(cid==j,3),'.')
            scatter3(r_brest(cid==j,1),r_brest(cid==j,2),r_brest(cid==j,3),'.','MarkerEdgeAlpha',0.65)
        end
    %     plot3(r_brest(cid==-1,1),r_brest(cid==-1,2),r_brest(cid==-1,3),'k.');
        scatter3(r_brest(cid==-1,1),r_brest(cid==-1,2),r_brest(cid==-1,3),5,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.03);
        hold off;
        title('knn-search')

        % Save figure
        if exist(fullfile(sessions(i).folder,'Figs_Tatsumi'),'dir')
            saveas(gcf,fullfile(sessions(i).folder,'Figs_Tatsumi','resting_clustering_knn_all.png'))
        else
            mkdir(fullfile(sessions(i).folder,'Figs_Tatsumi'))
            saveas(gcf,fullfile(sessions(i).folder,'Figs_Tatsumi','resting_clustering_knn_all.png'))
        end
        
        close gcf
        
        %=== SAVE: cluster index of each position
        save(fullfile(sessions(i).folder,'bpos.mat'),'bpos')
    
    end
    
    %=== FIGURE: clustering of resting positions for each bat
    figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0.2 1 0.38]);
    tl = tiledlayout(1, n_tags, 'TileSpacing', 'tight'); title(tl,sprintf('Clustering by DBSCAN & knn-search (day %d)',i))
    for j=1:n_tags
        nexttile(); hold on; box on; ax = gca; ax.BoxStyle = 'full';
        xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:)); xlabel("x"); ylabel("y"); zlabel("z"); view([135 50]);
        title(sprintf('%s (clustered=%.1f%%)',bat_nms(j,:),100*sum(~ismember(bpos(:,j),[-1 0]))/sum(~ismember(bpos(:,j),[0]))),'FontSize',12,'Color',bat_clr(j,:))
        for k=1:length(unique(bpos))-2 % exclude outliers (-1) and flying (0)
            scatter3(r(bpos(:,j)==k,1,j),r(bpos(:,j)==k,2,j),r(bpos(:,j)==k,3,j),'.','MarkerEdgeAlpha',0.65)
        end
        scatter3(r(bpos(:,j)==-1,1,j),r(bpos(:,j)==-1,2,j),r(bpos(:,j)==-1,3,j),5,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.03);
        hold off;
    end
    % Save figure
    if exist(fullfile(sessions(i).folder,'Figs_Tatsumi'),'dir')
        saveas(gcf,fullfile(sessions(i).folder,'Figs_Tatsumi','resting_cluster_knn.png'))
    else
        mkdir(fullfile(sessions(i).folder,'Figs_Tatsumi'))
        saveas(gcf,fullfile(sessions(i).folder,'Figs_Tatsumi','resting_cluster_knn.png'))
    end
    close gcf
    
    pause(0.3)
end

end