function PlotClusteredTrajectory(TrjCluster,mtData,FigDir,varargin)
% PlotClusterTrajectory(hc) plot the trajectories for each cluster identified by hierarchical clustering.
% Input
% hc: the output of flightClustering.m
% bat idx: a vector of bat index. The plots are generated only for the
% clusters for the bat specified by the bat idx vector.

p = inputParser;
addRequired(p, 'TrjCluster');
addRequired(p, 'FigDir')
addOptional(p, 'bat_idx', []);


parse(p, TrjCluster, FigDir, varargin{:});
bat_idx = p.Results.bat_idx; % index of bats
bnames = mtData.bname;

% Generate figures only for the specified bats
if isempty(bat_idx)
    bat_idx = 1:length(TrjCluster);
end

% Make directory for saving figures

%===FIGURE: visualize clustered trajectories for each bat
for i = 1:length(bat_idx)
    bb = bat_idx(i); % bat index
    nClu = TrjCluster(bb).n_clu;

    % plot all trajectories for each cluster
    fig = figure('Visible','off');
    set(fig, 'units','normalized','Position',[0 0 1 1]);
    tl = tiledlayout(3,6,"TileSpacing",'compact');
    title(tl,sprintf('Bat %s',bnames{bb}))

    for clu = 1:nClu
        clu_id = TrjCluster(bb).clu_list(clu); % index of cluster
        trj_count = TrjCluster(bb).clu_sz(clu); % count of trajectories in the current cluster
        clu_locs = TrjCluster(bb).clu_idx == clu_id; % locs of the current cluster
        r_flight = TrjCluster(bb).r(clu_locs); % trajectory of the current cluster
    
        nexttile

        % figure('Visible','off')
        % set(gcf, 'units', 'normalized', 'outerposition', [0.25 0.25 0.5 0.5]);
        hold on
        for trj = 1:trj_count
            x = r_flight{trj}(:,1);
            y = r_flight{trj}(:,2);
            z = r_flight{trj}(:,3);
            n = length(x);
            GradientColor = [uint8(jet(n)*255) uint8(ones(n,1))].';
            p = plot3(x,y,z);
            drawnow
            set(p.Edge, 'ColorBinding','interpolated', 'ColorData',GradientColor)
        end
        hold off
        title(sprintf('Cluster %d: %d flights',clu_id,trj_count))
        xlabel('x')
        ylabel('y')
        zlabel('z')
        xlim([-3,3])
        ylim([-3,3])
        zlim([0,2.5])

        view([0,90])
    end
    if ~isempty(FigDir)
        if isempty(dir(FigDir))
            mkdir(FigDir)
        end
        saveas(fig,fullfile(FigDir,sprintf('TrajectoryCluster_bat%s.jpg',bnames{bb})))
    end
    close(fig)
end
end