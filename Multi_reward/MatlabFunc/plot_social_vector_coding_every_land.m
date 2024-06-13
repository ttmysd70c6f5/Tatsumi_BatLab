function plot_social_vector_coding_every_land(BhvData,TrjCluster,spikeData,mtData,saveDir)
n_bat = BhvData.n_tags;
n_unit = spikeData.numGoodUnits;
unit_ids = spikeData.good_units.cluster_id;
regions = spikeData.good_units.region;

%%
load(fullfile(saveDir,'vector_coding_every_pos.mat'),"Xedges",'Yedges','social_vector_map','bb')
%%% FIGURE: agent-vector coding for each landing position
bnames = mtData.bname;
for bbb = [1:bb-1,bb+1:n_bat]
    clu_list = social_vector_map(bbb).clu_list;
    rate_map = social_vector_map(bbb).rate_map;
    occupancy_map_flt = social_vector_map(bbb).occupancy_map_flt;
        
    for unit = 1:n_unit
        unit_id = unit_ids(unit);
        region = regions{unit};

        % fig = figure;
        fig = figure('Visible','off');
        set(fig, 'units','normalized','Position',[0 0.3 1 0.3]);
        tl = tiledlayout(1,length(clu_list),"TileSpacing",'compact');
        title(tl,sprintf('%s-%s, Unit%d, %s',bnames{bb},bnames{bbb},unit_id,region),'Interpreter','none')

        for clu = 1:length(clu_list)
            clu_id = clu_list(clu);
            r_land = TrjCluster(bbb).r_land_mean(TrjCluster(bbb).clu_list == clu_id,1:2); % target bat position



            % fig = figure('visible','off');
            nexttile
            x = Xedges(1:end-1) + diff(Xedges)/2;
            y = Yedges(1:end-1) + diff(Yedges)/2;
            peak_fr = max(rate_map(:,:,clu,unit),[],'all');
            hold on 
            h = imagesc([x(1) x(end)],[y(1) y(end)],rate_map(:,:,clu,unit));
            set(gca,'YDir','normal','FontSize',16)
            set(h,'AlphaData',~isnan(occupancy_map_flt(:,:,clu)))

            scatter(r_land(1),r_land(2),50,'r','o','filled')
            hold off

            % ylabel(sprintf('%s_%d',bnames{bb},unit_id))
            % ylabel('Y (m)')
            title(sprintf('Clu %d, %.1f Hz',clu_id,peak_fr),'Interpreter','none')
            set(gca,'XTick',[],'YTick',[])
            colormap parula
            colorbar

            % % Save
            % figDir = fullfile(saveDir,bnames{bbb},sprintf('unit%d',unit_id));
            % if isempty(dir(figDir)); mkdir(figDir); end
            % saveas(fig,fullfile(figDir,sprintf('clu%d.png',clu_id)))
            % close(fig)
        end
        % Save
        figDir = fullfile(saveDir,bnames{bbb});
        if isempty(dir(figDir)); mkdir(figDir); end
        saveas(fig,fullfile(figDir,sprintf('unit%d.png',unit_id)))
        close(fig)
    end
end

%%% FIGURE: allocentric spatial occupancy
bnames = mtData.bname;
for bbb = [1:bb-1,bb+1:n_bat]
    clu_list = social_vector_map(bbb).clu_list;
    occupancy_map_flt = social_vector_map(bbb).occupancy_map_flt;
        
    fig = figure;
    % fig = figure('Visible','off');
    set(fig, 'units','normalized','Position',[0 0.3 1 0.3]);
    tl = tiledlayout(1,length(clu_list),"TileSpacing",'compact');
    title(tl,sprintf('%s-%s: spatial occupancy (sec)',bnames{bb},bnames{bbb}),'Interpreter','none')

    for clu = 1:length(clu_list)
        clu_id = clu_list(clu);

        % fig = figure('visible','off');
        nexttile
        x = Xedges(1:end-1) + diff(Xedges)/2;
        y = Yedges(1:end-1) + diff(Yedges)/2;
        h = imagesc([x(1) x(end)],[y(1) y(end)],occupancy_map_flt(:,:,clu));
        set(gca,'YDir','normal','FontSize',16)
        set(h,'AlphaData',~isnan(occupancy_map_flt(:,:,clu)))
        xlabel('X (m)')
        ylabel('Y (m)')
        % title(sprintf('%s-%s: occupancy (s)',bnames{bb},bnames{bbb}),'Interpreter','none')
        title(sprintf('Clu %d',clu_id))
        colormap jet
        colorbar
    
        % % Save
        % figDir = fullfile(saveDir,'spatial_occupancy');
        % if isempty(dir(figDir)); mkdir(figDir); end
        % saveas(fig,fullfile(figDir,sprintf('%s_clu%d.png',bnames{bbb},clu_id)))
        % close(fig)
    end
    % Save
    figDir = fullfile(saveDir,bnames{bbb});
    if isempty(dir(figDir)); mkdir(figDir); end
    saveas(fig,fullfile(figDir,'Spatial_occupancy.png'))
    close(fig)
end