function flight_clus = FlightClus_AF_v3(r_vector,f_vector,Fs,varargin)
%Flight_clus_AF extracts flight clusters from trajectory data, by using
%hierarchical agglomerative clustering (or k-means) on 3D points along the trajectory.

%INPUTS:
% r_vector = time (rows) x 3D position (columns)
% f_vector = time (rows) x flight (column): 0s (not flying) or 1s (flying)
% varargin =

%OUTPUTS:
% flight_clus STRUCTURE
% flight_clus.id:           id of the cluster (cluster 1 pools all unclustered)                
% flight_clus.strt_frame:   takeoff frame
% flight_clus.stop_frame:   landing frame
% flight_clus.pos:          3 x samples x N matrix with all the coordinates
% flight_clus.vel:          1 x samples x N matrix with all the velocities
% flight_clus.Fs:           Sampling Frequency
% flight_clus.N:            Total number of flights
% flight_clus.length:       Lenght of the flight
% flight_clus.dur:          Duration of the flight
% flight_clus.sorted_ids:   Ids, sorted according to takeoff sample

% COMMENTS:
%210912: the function is pretty sensitive to the Flight Segmentation
%210912: 6 pts,euclidean distance; 10 pts, PCA

% USAGE EXAMPLES
% FlightClus_AF_v2(squeeze(r(:,:,i)),bflying(:,i),'Alpha',alpha_clus,'Frechet',1,'Points',10);
% FlightClus_AF_v2(squeeze(r(:,:,i)),bflying(:,i),'Alpha',alpha_clus,'Frechet',0,'Points',6);


%% Default Parameters and overrides
%x2=2.9; x1=-2.9; y2=2.6;  y1=-2.6;  z1=0; z2=2.3;  %flight volume coordinates
%x2=5.5; x1=-4.6; y2=1.8;  y1=-2.3;  z1=0; z2=2.8;  %flight volume coordinates
x2 = max(r_vector(:,1))*1.1;    x1 = min(r_vector(:,1))*1.1;
y2 = max(r_vector(:,2))*1.1;    y1 = min(r_vector(:,2))*1.1;
z2 = max(r_vector(:,3))*1.1;    z1 = min(r_vector(:,3))*1.1;
ds_clus = 10;                                       %number of 3D-points/flight for clustering (depending of the clustering algorithm, usually >6 points is enough)
pca_features = 0;                                   %if using PCA
Frechet_clus = 0;                                   %if using Frechet-distance for clustering
k_means = 0;                                        %if using k-means
alpha = 3.5;                                        %main clustering parameter
reassign = 1;                                       %re-order clusters
N_min = 3;                                          %min number of flights for being a cluster
Frechet = 0;                                        %if performing Frechet-based calculations
calc_length = 1;                                    %length and duration calculation (takes time)
ignore_z = 0;                                       %if projecting everithing on the xy plane

%User input overrides default params
if nargin > 1
    nparams=length(varargin);
    for i=1:2:nparams
        switch (varargin{i})
            case 'Fs'
                Fs=varargin{i+1};
            case 'Points'
                ds_clus = varargin{i+1};
            case 'N_min'
                N_min = varargin{i+1};
            case 'Alpha'
                alpha = varargin{i+1};
            case 'Frechet'
                Frechet_clus = varargin{i+1};
            case 'PCA'
                pca_features = varargin{i+1};
        end
    end
end

%Sanity check: a flight is not cut at the beginning or end of the session
if f_vector(1) == 1 || f_vector(end) == 1
    warning('INCOMPLETE FLIGHTS AT START/STOP OF THE SESSION');
end

%Extract flight number, start and stop and calculate velocity
num_flights = nnz(diff(f_vector)>0);
f_start = find(diff(f_vector)>0)+1;
f_stop = find(diff(f_vector)<0);
v = diff(r_vector,1,1).*Fs; v=[zeros(1,3); v];   v_abs = vecnorm(v,2,2);


%% Define features

%Cut out flights, downsample to ds_clus positions per flight
all_flights = NaN(3,max(f_stop-f_start)+1,num_flights);          %3D matrix with all flights (position)
all_flights_vel = NaN(1,max(f_stop-f_start)+1,num_flights);      %3D matrix with all flights (velocity)
all_flights_ds = NaN(3,ds_clus,num_flights);                     %3D matrix with all flights(downsampled position)

for nf = 1 : num_flights
    all_flights(:,1:(f_stop(nf)-f_start(nf)+1),nf) = r_vector(f_start(nf):f_stop(nf),:)';
    all_flights_vel(1,1:(f_stop(nf)-f_start(nf)+1),nf) = v_abs(f_start(nf):f_stop(nf),:)';
    all_flights_ds(:,:,nf) = interp1(linspace(1,100,f_stop(nf)-f_start(nf)+1)',r_vector(f_start(nf):f_stop(nf),:),linspace(1,100,ds_clus)','linear')';
    
    %     %Uncomment if you want to see how the downsampled flights look like
    %     h = figure();   set(gcf, 'units','normalized','outerposition',[0.3 0.5 0.4 0.4]);
    %     subplot(3,5,[1,2,6,7,11,12]);
    %     plot3(all_flights(1,:,nf),all_flights(2,:,nf),all_flights(3,:,nf),'Color','b');   hold on;
    %     plot3(all_flights_ds(1,:,nf),all_flights_ds(2,:,nf),all_flights_ds(3,:,nf),'.','Color','r','MarkerSize',10);  hold off;
    %     subplot(3,5,[3:5]);     plot(linspace(1,100,f_stop(nf)-f_start(nf)+1)',all_flights(1,~isnan(all_flights(1,:,nf)),nf),linspace(1,100,ds_clus)',all_flights_ds(1,:,nf),'.');
    %     subplot(3,5,[8:10]);    plot(linspace(1,100,f_stop(nf)-f_start(nf)+1)',all_flights(2,~isnan(all_flights(2,:,nf)),nf),linspace(1,100,ds_clus)',all_flights_ds(2,:,nf),'.');
    %     subplot(3,5,[13:15]);   plot(linspace(1,100,f_stop(nf)-f_start(nf)+1)',all_flights(3,~isnan(all_flights(3,:,nf)),nf),linspace(1,100,ds_clus)',all_flights_ds(3,:,nf),'.');
    %     waitfor(h);     %w = waitforbuttonpress;
end

%Define X matrix of features for clustering (downsampled coordinates, stacked together)
if ~Frechet_clus
    X = reshape(permute(all_flights_ds,[2 1 3]),[],num_flights)';   %so then X = #flights x #features
else
    X = all_flights_ds;
end

if ignore_z,X(3,:,:) = 0;end

%Define weights to increase the importance of takeoff and landing spots
%e_weights = 1-0.4*gausswin(ds_clus);
%e_weights = repmat(e_weights,3,1);
e_weights = repmat(ones(ds_clus,1),3,1);

%If dimensionality reduction is needed
if pca_features
    [~,score,~,~,explained,~] = pca(X);     X = score(:,1:4);
end

%% Perform clustering
Y = [];
if k_means
    n_clusters = 15;    idx = kmeans(X,n_clusters);
else
    dist = alpha;                               %linkage distance
    if ~Frechet_clus
        Y = pdist(X,'seuclidean',e_weights);
        Y = Y./ds_clus;
    else
        pairs = nchoosek(1:num_flights,2);
        for n = 1:length(pairs)
            Y(n) = DiscreteFrechetDist(squeeze(X(:,:,pairs(n,1))),squeeze(X(:,:,pairs(n,2))));
        end
    end
    Z = linkage(Y,'single');
    idx = cluster(Z,'Cutoff',dist,'Criterion','distance');
end

%% Create structure with features of the clusters and rearrange

flight.strt_frame = ceil(f_start)';
flight.stop_frame = ceil(f_stop)';
flight.pos = all_flights;
flight.vel = all_flights_vel;
flight.id = idx;
flight.Fs = Fs;

%Sort structure according to cluster id
[flight_sorted.id,I] = sort(flight.id);
flight_sorted.strt_frame = flight.strt_frame(I);
flight_sorted.stop_frame = flight.stop_frame(I);
flight_sorted.pos = flight.pos(:,:,I);
flight_sorted.vel = flight.vel(:,:,I);
flight_sorted.Fs = flight.Fs;
flight_sorted.N = size(flight_sorted.id,1);

%Assign isolated clusters to cluster #flights+1
[Ns,b] = histc(flight_sorted.id,unique(flight_sorted.id));
flight_sorted.id(Ns(b)<N_min) = size(all_flights,3)+1;
id_surv_clusters = unique(flight_sorted.id);
n_surv_clusters = size(id_surv_clusters,1);

%Create final structure flight.clus after re-assignment
%unclustered flights are in the last cluster
flight_clus.id = flight_sorted.id;
flight_clus.strt_frame = flight_sorted.strt_frame;
flight_clus.stop_frame = flight_sorted.stop_frame;
flight_clus.pos = flight_sorted.pos;
flight_clus.vel = flight_sorted.vel;
flight_clus.Fs = flight_sorted.Fs;
flight_clus.N = flight_sorted.N;
for jj=1:n_surv_clusters
    flight_clus.id(flight_sorted.id == id_surv_clusters(jj)) = jj;
end
id_surv_clusters = unique(flight_clus.id);

%Re-assign id for convenience, with cluster 1 of non-clustered flights and
%then clusters reordered in descending order of crowdedness
if reassign
    new_ord = [];   [~,new_ord] = sort(histc(flight_clus.id,id_surv_clusters(1:end-1)),'descend');
    new_ord = [new_ord; id_surv_clusters(end)];
    new_ord = circshift(new_ord,1);
    reassign_matrix =(flight_clus.id == new_ord');
    for jj=1:n_surv_clusters
        flight_clus.id(reassign_matrix(:,jj)) = jj;
    end
end

%Calculate trajectory length, linear parametrization, duration in s and interflight (take-off to take-off)
for ii = 1:flight_sorted.N
    if calc_length
%         % More precise method, based on spline interpolation, need to get rid of replicates
%         trajectory = unique(round(flight_clus.pos(:,~isnan(flight_clus.pos(1,:,ii)),ii)',3),'rows','stable');
%         [flight_clus.length(ii), seglen] = arclength(trajectory(:,1),trajectory(:,2),trajectory(:,3),'s');
        % Less precise method, based on linear interpolation, faster but biased
        trajectory = flight_clus.pos(:,~isnan(flight_clus.pos(1,:,ii)),ii)';
        [flight_clus.length(ii), seglen] = arclength(trajectory(:,1),trajectory(:,2),trajectory(:,3),'l');
        flight_clus.lin_tr{1,ii} = [0; cumsum(seglen)]/flight_clus.length(ii);
        flight_clus.dur(ii) = (flight_clus.stop_frame(ii)-flight_clus.strt_frame(ii))./flight_clus.Fs;
    else
        flight_clus.length(ii)= 0;
        flight_clus.dur(ii) = 0;
    end
end

% Provide Cluster ids sorted from first to last flight (OBSOLETE)
[~,I] = sort(flight_clus.strt_frame);
sorted_ids = flight_clus.id(I);
flight_clus.sorted_ids = sorted_ids;

% Sort Structure based on Take-off sample, Fix Transpositions and ds_pos
[~,I] = sort(flight_clus.strt_frame);
flight_clus.id = flight_clus.id(I,:)';
flight_clus.strt_frame = flight_clus.strt_frame(:,I);
flight_clus.stop_frame = flight_clus.stop_frame(:,I);
flight_clus.pos = flight_clus.pos(:,:,I);
flight_clus.vel = flight_clus.vel(:,:,I);
flight_clus.length = flight_clus.length(:,I);
flight_clus.lin_tr = {flight_clus.lin_tr{:,I'}};
flight_clus.dur = flight_clus.dur(:,I);
flight_clus.sorted_ids = flight_clus.sorted_ids';
flight_clus.ds_pos = all_flights_ds;

% Calculate correlation and Frechet distance with mean of the cluster
flight_clus.corr = zeros(1,num_flights);    flight_clus.dist = zeros(1,num_flights);
for jj=1:n_surv_clusters
    mean_path = mean(flight_clus.ds_pos(:,:,find(flight_clus.id == jj)),3);
    for i = find(flight_clus.id == jj)
        flight_clus.corr(:,i) = mean(spdiags(corr(flight_clus.ds_pos(:,:,i)',mean_path'),0));
        [flight_clus.dist(:,i),~] = DiscreteFrechetDist(flight_clus.ds_pos(:,:,i)',mean_path');
    end
end

% Be careful: correlation does not capture offsets or variations in the
% flight 'amplitude' (think at a sine wave)

% %=== Visualize flights in cluster jj
% figure('units','normalized','outerposition',[0.6 0.2 0.2 0.5]);
% for jj=5
%     id = find(flight_clus.id==jj);
%     plot3(r_vector(:,1),r_vector(:,2),r_vector(:,3),':','Color',[0.8 0.8 0.8],'MarkerSize',0.001);
%     xlim([x1 x2]); ylim([y1 y2]);   zlim([z1 z2]);  view(0,90);
%     xlabel('x');    ylabel('y');
%     hold on;
%     mean_path = mean(flight_clus.ds_pos(:,:,find(flight_clus.id == jj)),3);
%     plot3(mean_path(1,:),mean_path(2,:),mean_path(3,:),'-','LineWidth',3,'Color', 'k');
%     for ii=1:size(id,2)
%         title(['Cluster' num2str(jj) ' (' num2str(size(id,2)) ' flights),' num2str(size(id,2)/num_flights*100,2) '%'])
%         plot3(flight_clus.pos(1,:,id(ii)),flight_clus.pos(2,:,id(ii)),flight_clus.pos(3,:,id(ii)),'-','LineWidth',1.5,'Color', col(jj,:));
%         text(flight_clus.ds_pos(1,round(ds_clus/2),id(ii)),flight_clus.ds_pos(2,round(ds_clus/2),id(ii)),flight_clus.ds_pos(3,round(ds_clus/2),id(ii)),num2str(flight_clus.corr(:,id(ii)),3));
%     end
%     hold off;
% end

%% !!!!--------> PLOTTING <--------!!!! 

%=== Cluster Colors
col = hsv(n_surv_clusters);

%=== Embed in 2-D space
f_embd = ones(flight_clus.N,2);
%f_embd = mdscale(Y,2,'Criterion','sammon'); 

%=== Plot distances between flights and dendrogram, plot embedding with clusters
figure('units','normalized','outerposition',[0.5 0.1 0.15 0.8]);
subplot(411); edges_f = 10.^linspace(log10(min(Y)),log10(max(Y)),50);  histogram(Y,edges_f,'EdgeColor','none');    set(gca,'XScale','log');    set(gca,'YScale','log'); xlabel('Inter-flight distances (m)');
subplot(412); histogram(Y,100,'EdgeColor','none'); xlabel('Inter-flight distances (m)');
subplot(413); dendrogram(Z,0);  hold on;    refline(0,dist);    hold off;
subplot(414); gscatter(f_embd(:,1),f_embd(:,2),flight_clus.id,[],[],10,'off');  axis equal;
        
% %=== Visualize clusters in PCA space
% figure;
% if ~pca_features && ~Frechet_clus
%     [~,score,~,~,explained,~] = pca(X);
% elseif ~pca_features && Frechet_clus
%     [~,score,~,~,explained,~] = pca(reshape(permute(all_flights_ds,[2 1 3]),[],num_flights)');
% end
% hold on;
% for jj=1:n_surv_clusters
%     id = find(sorted_ids==jj);
%     if jj == 1
%         sz = 10;
%     else
%         sz = 36;
%     end
%     scatter3(score(id,1),score(id,2),score(id,3),sz,'o','MarkerFaceColor',col(jj,:),'MarkerEdgeColor',[1 1 1]);
%     xlabel('PC1');   ylabel('PC2');   zlabel('PC3');
%     var_expl = cumsum(explained);   title([num2str(var_expl(3),2) '% Variance explained by the first 3 PCs'])
% end
% hold off;

if 1
    
    %=== Visualize flights in cluster 1
    figure('units','normalized','outerposition',[0.6 0.2 0.2 0.5]);
    for jj=1
        id = find(flight_clus.id==jj);
        plot3(r_vector(:,1),r_vector(:,2),r_vector(:,3),':','Color',[0.8 0.8 0.8],'MarkerSize',0.001);
        xlim([x1 x2]); ylim([y1 y2]);   zlim([z1 z2]);  view(0,90);
        xlabel('x');    ylabel('y');
        hold on;
        for ii=1:size(id,2)
            title(['Cluster' num2str(jj) ' (' num2str(size(id,2)) ' flights),' num2str(size(id,2)/num_flights*100,2) '%'])
            plot3(flight_clus.pos(1,:,id(ii)),flight_clus.pos(2,:,id(ii)),flight_clus.pos(3,:,id(ii)),'-','LineWidth',1.5,'Color', col(jj,:));
        end
        hold off;
    end
    
    %=== Plot All Clusters (1D)
    figure('units','normalized','outerposition',[0.25 0 0.3 1]);
    tiledlayout(n_surv_clusters,6,'TileSpacing','none','Padding','compact');
    titles = {'x','y','z','x DS','y DS','z DS'};
    for jj=1:n_surv_clusters
        id = find(flight_clus.id==jj);
        for c = 1:3
            nexttile;  hold on;
            if c==1, ylabel(['Cluster', num2str(jj)]);    end
            if jj==1,title(titles{1,c});end
            for ii=1:size(id,2)
                trace = flight_clus.pos(c,~isnan(flight_clus.pos(c,:,id(ii))),id(ii));
                plot(linspace(1,100,length(trace)),trace,'-','LineWidth',0.3,'Color', col(jj,:));
            end
            hold off;
        end
        for c = 1:3
            nexttile;   hold on;
            if jj==1,title(titles{1,c+3});end
            for ii=1:size(id,2)
                trace = all_flights_ds(c,~isnan(all_flights_ds(c,:,id(ii))),id(ii));
                plot(linspace(1,100,length(trace)),trace,'-','LineWidth',0.3,'Color', col(jj,:));
            end
            hold off;
        end
    end
    
    %=== Plot All Clusters (3D)
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(3,n_surv_clusters,'TileSpacing','tight','Padding','compact');
    for jj=1:n_surv_clusters
        id = find(flight_clus.id==jj);
        
        nexttile(0*n_surv_clusters+jj);
        plot3(r_vector(:,1),r_vector(:,2),r_vector(:,3),':','Color',[0.8 0.8 0.8],'MarkerSize',0.001);
        hold on;
        avg_take_off = [];
        for ii=1:size(id,2)
            title(['Cluster' num2str(jj) ' (' num2str(size(id,2)) ' flights),' num2str(size(id,2)/num_flights*100,2) '%'])
            plot3(flight_clus.pos(1,:,id(ii)),flight_clus.pos(2,:,id(ii)),flight_clus.pos(3,:,id(ii)),'-','LineWidth',1.5,'Color', col(jj,:));
            %text(flight_clus.ds_pos(1,round(ds_clus/2),id(ii)),flight_clus.ds_pos(2,round(ds_clus/2),id(ii)),num2str(flight_clus.dist(:,id(ii)),3));
            avg_take_off = [avg_take_off flight_clus.pos(:,1,id(ii))];
        end
        take_off = mean(avg_take_off,2);    textscatter(take_off(1),take_off(2),"Take-off");        hold off;
        axis equal; xlim([x1 x2]); ylim([y1 y2]);   zlim([z1 z2]);  view(0,90); xlabel('x');    ylabel('y');    
        
        nexttile(1*n_surv_clusters+jj);
        plot(r_vector(:,1),r_vector(:,3),':','Color',[0.8 0.8 0.8],'MarkerSize',0.001);
        hold on;
        avg_take_off = [];
        for ii=1:size(id,2)
            plot(flight_clus.pos(1,:,id(ii)),flight_clus.pos(3,:,id(ii)),'-','LineWidth',1.5,'Color', col(jj,:));
            %text(flight_clus.ds_pos(1,round(ds_clus/2),id(ii)),flight_clus.ds_pos(3,round(ds_clus/2),id(ii)),num2str(flight_clus.dist(:,id(ii)),3));
            avg_take_off = [avg_take_off flight_clus.pos(:,1,id(ii))];
        end
        take_off = mean(avg_take_off,2);    textscatter(take_off(1),take_off(3),"Take-off");     hold off;
        axis equal; xlim([x1 x2]); ylim([z1 z2]);       xlabel('x');    ylabel('y');    
        
        nexttile(2*n_surv_clusters+jj);
        histogram(flight_clus.dur(id),'EdgeColor','none'); xlabel('Duration(s)');  ylabel('Counts');
    end
end

%=== Plot flights stats
% figure;
% subplot(1,3,1);     histogram(flight_clus.length);   xlabel('Flight length (m)');        ylabel('Counts');
% subplot(1,3,2);     histogram(flight_clus.dur);      xlabel('Fligth duration (s)');      ylabel('Counts');

%Calculate entropy of the session
%[Pr,b] = histc(flight.id,unique(flight.id));    Pr = Pr./flight_clus.N;     Entropy = -sum(Pr.*log2(Pr));

%% Calculate the Frechet distance between flights
if Frechet
    %Calculate distance between all the flights
    flight_pairs = nchoosek(1:flight_clus.N,2);
    F_dist = zeros(size(flight_pairs,1),1);
    for i = 1:size(flight_pairs,1)
        P = flight_clus.pos(:,:,flight_pairs(i,1));     P = P(:,~isnan(P(1,:)))';   P = downsample(P,3);
        Q = flight_clus.pos(:,:,flight_pairs(i,2));     Q = Q(:,~isnan(Q(1,:)))';   Q = downsample(Q,3);
        [F_dist(i),~] = DiscreteFrechetDist(P,Q);
    end
    
    %Intra-cluster Frechet distance
    [Nf,~] = histc(flight_clus.id,unique(flight_clus.id));
    F_dist_intra = nan(max(Nf)*(max(Nf)-1)/2,n_surv_clusters);
    if n_surv_clusters>1
        for jj = 1:n_surv_clusters
            disp(jj);
            flight_pairs = nchoosek(find(flight_clus.id==jj),2);
            for i = 1:size(flight_pairs,1)
                P = flight_clus.pos(:,:,flight_pairs(i,1));     P = P(:,~isnan(P(1,:)))';   P = downsample(P,3);
                Q = flight_clus.pos(:,:,flight_pairs(i,2));     Q = Q(:,~isnan(Q(1,:)))';   Q = downsample(Q,3);
                [F_dist_intra(i,jj),~] = DiscreteFrechetDist(P,Q);
            end
        end
    end
    
    %Extra-cluster Frechet distance
    if n_surv_clusters>1
        flight_pairs = nchoosek(1:flight_clus.N,2);
        idx_del = false(size(flight_pairs,1),1);
        for jj = 2:n_surv_clusters
            to_delete = nchoosek(find(flight_clus.id==jj),2);
            idx_del = idx_del | ismember(flight_pairs,to_delete,'rows') | ismember(flight_pairs,to_delete(:,[2 1]),'rows');
        end
        F_dist_extra = F_dist(~idx_del);
    else
        F_dist_extra = F_dist;
    end
    
    %Stats for different clusters (cv duration, length, Frechet and interflight interval for same cluster)
    Clus = [];
    for jj = 1:n_surv_clusters
        Clus.cv_dur(jj) = std(flight_clus.dur(ismember(flight_clus.id,jj)))/mean(flight_clus.dur(ismember(flight_clus.id,jj)));
        Clus.cv_len(jj) = std(flight_clus.length(ismember(flight_clus.id,jj)))/mean(flight_clus.length(ismember(flight_clus.id,jj)));
        Clus.cv_fre(jj) = mean(F_dist_intra(:,jj),'omitnan')/mean(flight_clus.length(ismember(flight_clus.id,jj)));
        Clus.if_int(jj) = median(diff(flight_clus.strt_frame(ismember(flight_clus.id,jj)))./flight_clus.Fs);
    end
    
    figure('units','normalized','outerposition',[0.1 0.25 0.8 0.5]);
    ax_1 = subplot(1,3,1);
    hfl = histogram(F_dist,'facealpha',.5,'edgecolor','none');
    xlabel('Frechet distance (m)'); ylabel('Counts'); title('All flights');
    
    ax_2 = subplot(1,3,2);
    for jj = 1:n_surv_clusters
        histogram(F_dist_intra(:,jj),hfl.BinEdges,'facecolor',col(jj,:),'facealpha',.5,'edgecolor','none'); hold on;
    end
    hold off;   xlabel('Frechet distance (m)'); title('Intra clusters 1,2,...');
    
    ax_3 = subplot(1,3,3);
    histogram(F_dist_extra,hfl.BinEdges,'facecolor','b','facealpha',.5,'edgecolor','none');   hold on;
    histogram(F_dist_intra(:,2:end),hfl.BinEdges,'facecolor','k','facealpha',.5,'edgecolor','none');   hold off;
    xlabel('Frechet distance (m)');            title('Intra VS Extra');
    
    linkaxes([ax_1,ax_2,ax_3],'y');
    %set([ax_1,ax_2,ax_3], 'YScale', 'log');
    
end
    
end