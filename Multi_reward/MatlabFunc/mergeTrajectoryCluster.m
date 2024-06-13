function TrjCluster = mergeTrajectoryCluster(TrjCluster,FigDir,bid,old_cluster_1,old_cluster_2,varargin)

% current cluster
cluster = TrjCluster(bid).cluster;
cluster_count = TrjCluster(bid).count;
cluster_idx = TrjCluster(bid).cluster_idx;
cluster_idx_max = max(cluster_idx);
num_cluster = TrjCluster(bid).num_cluster;
r_mean_toff = TrjCluster(bid).r_mean_toff;
r_mean_tland = TrjCluster(bid).r_mean_land;
r_mean_interp = TrjCluster(bid).r_interp_mean;

% update cluster
old_id_1 = find(cluster_idx == old_cluster_1);
old_id_2 = find(cluster_idx == old_cluster_2);
if nargin == 5
    new_cluster_id = cluster_idx_max + 1; % newly assigned cluster id
elseif nargin == 6
    new_cluster_id = varargin{:};
else
    disp("Incorrect number of input arguments")
end

cluster(cluster==old_cluster_1) = new_cluster_id;
cluster(cluster==old_cluster_2) = new_cluster_id;
cluster_count(end+1) = cluster_count(old_id_1) + cluster_count(old_id_2);
cluster_count(old_id_1) = 0;
cluster_count(old_id_2) = 0;
cluster_idx(end+1) = new_cluster_id;
num_cluster = num_cluster - 1;
r_mean_toff(end+1,:) = mean([r_mean_toff(old_id_1,:); r_mean_toff(old_id_2,:)],1);
r_mean_tland(end+1,:) = mean([r_mean_tland(old_id_1,:); r_mean_tland(old_id_2,:)],1);
r_mean_interp(end+1,:) = mean([r_mean_interp(old_id_1,:); r_mean_interp(old_id_2,:)],1);

% plot new trajectory
figure('Visible','off')
% set(gcf, 'units', 'normalized', 'outerposition', [0.25 0.25 0.5 0.5]);
hold on

r_flight = TrjCluster(bid).r_flight(cluster==new_cluster_id);
n_trj = cluster_count(cluster_idx==new_cluster_id);
for trj = 1:n_trj
    x = r_flight{trj}(:,1,bid);
    y = r_flight{trj}(:,2,bid);
    z = r_flight{trj}(:,3,bid);
    n = length(x);
    GradientColor = [uint8(jet(n)*255) uint8(ones(n,1))].';
    p = plot3(x,y,z);
    drawnow
    set(p.Edge, 'ColorBinding','interpolated', 'ColorData',GradientColor)
end
hold off
title(sprintf('Bat %d, cluster %d',bid, new_cluster_id))
xlabel('x')
ylabel('y')
zlabel('z')
xlim([-3,3])
ylim([-3,3])
zlim([0,2.5])

view([0,90])

if ~isempty(FigDir)
    if isempty(dir(FigDir))
        mkdir(FigDir)
    end
    saveas(gcf,fullfile(FigDir,sprintf('TrajectoryCluster_bat%d_cluster%d.jpg',bid,new_cluster_id)))
end
close gcf


% save to struct
TrjCluster(bid).cluster = cluster;
TrjCluster(bid).count = cluster_count;
TrjCluster(bid).cluster_idx = cluster_idx;
TrjCluster(bid).num_cluster = num_cluster - 1;
TrjCluster(bid).r_mean_toff = r_mean_toff;
TrjCluster(bid).r_mean_land = r_mean_tland;
TrjCluster(bid).r_interp_mean = r_mean_interp;

fprintf('Bat %d: Cluster %d and %d was merged into the cluster %d\n',bid,old_cluster_1,old_cluster_2,new_cluster_id)
end