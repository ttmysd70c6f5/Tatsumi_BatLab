function PlotClusteredTrajectory(TrjCluster,varargin)
% PlotClusterTrajectory(hc) plot the trajectories for each cluster identified by hierarchical clustering.
% Input
% hc: the output of flightClustering.m
% bat idx: a vector of bat index. The plots are generated only for the
% clusters for the bat specified by the bat idx vector.

p = inputParser;
addRequired(p, 'TrjCluster');
addOptional(p, 'bat_idx', []); % all = all landings, fd = only feeder landings
addOptional(p, 'FigDir', [])

parse(p, TrjCluster, varargin{:});
bat_idx = p.Results.bat_idx; % index of bats
FigDir = p.Results.FigDir; % the path to save figures


% Generate figures only for the specified bats
if isempty(bat_idx)
    bat_idx = 1:length(TrjCluster);
end

% Make directory for saving figures

%===FIGURE: visualize clustered trajectories for each bat
for i = 1:length(bat_idx)
    bb = bat_idx(i); % bat index
    nClu = TrjCluster(bb).num_cluster;
    for clu = 1:nClu
        cluster_idx = TrjCluster(bb).cluster_idx(clu); % index of cluster
        trj_count = TrjCluster(bb).count(clu); % count of trajectories in the current cluster
        clu_locs = TrjCluster(bb).cluster == cluster_idx; % locs of the current cluster
        r_flight = TrjCluster(bb).r_flight(clu_locs); % trajectory of the current cluster
    
        % plot all trajectories for each cluster
        figure('Visible','off')
        % set(gcf, 'units', 'normalized', 'outerposition', [0.25 0.25 0.5 0.5]);
        hold on
        for trj = 1:trj_count
            x = r_flight{trj}(:,1,bb);
            y = r_flight{trj}(:,2,bb);
            z = r_flight{trj}(:,3,bb);
            n = length(x);
            GradientColor = [uint8(jet(n)*255) uint8(ones(n,1))].';
            p = plot3(x,y,z);
            drawnow
            set(p.Edge, 'ColorBinding','interpolated', 'ColorData',GradientColor)
        end
        hold off
        title(sprintf('Bat %d, cluster %d',bb, cluster_idx))
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
            saveas(gcf,fullfile(FigDir,sprintf('TrajectoryCluster_bat%d_cluster%d.jpg',bb,cluster_idx)))
        end
        close gcf
    end
end
end