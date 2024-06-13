function plot_social_vector_coding(BhvData,spikeData,mtData,saveDir)
n_bat = BhvData.n_tags;
n_unit = spikeData.numGoodUnits;
unit_ids = spikeData.good_units.cluster_id;
regions = spikeData.good_units.region;

%%
load(fullfile(saveDir,'vector_coding.mat'),"Xedges",'Yedges','occupancy_map_flt','rate_map','bb')
%%% FIGURE: agent-vector coding
bnames = mtData.bname;
for bbb = [1:bb-1,bb+1:n_bat]
    for unit = 1:n_unit
        unit_id = unit_ids(unit);
        region = regions{unit};

        fig = figure('visible','off');
        x = Xedges(1:end-1) + diff(Xedges)/2;
        y = Yedges(1:end-1) + diff(Yedges)/2;
        peak_fr = max(rate_map(:,:,unit,bbb),[],'all');
        h = imagesc([x(1) x(end)],[y(1) y(end)],rate_map(:,:,unit,bbb));
        set(gca,'YDir','normal','FontSize',16)
        set(h,'AlphaData',~isnan(occupancy_map_flt(:,:,bbb)))
        xlabel('X (m)')
        ylabel('Y (m)')
        title(sprintf('%s-%s, unit %d, %s, %.1f Hz',bnames{bb},bnames{bbb},unit_id,region,peak_fr),'Interpreter','none')
        colormap jet
        colorbar

        % Save
        figDir = fullfile(saveDir,bnames{bbb});
        if isempty(dir(figDir)); mkdir(figDir); end
        saveas(fig,fullfile(figDir,sprintf('Unit%d.png',unit_id)))
        close(fig)
    end
end

%%% FIGURE: allocentric spatial occupancy
bnames = mtData.bname;
for bbb = [1:bb-1,bb+1:n_bat]
    fig = figure('visible','off');
    x = Xedges(1:end-1) + diff(Xedges)/2;
    y = Yedges(1:end-1) + diff(Yedges)/2;
    h = imagesc([x(1) x(end)],[y(1) y(end)],occupancy_map_flt(:,:,bbb));
    set(gca,'YDir','normal','FontSize',16)
    set(h,'AlphaData',~isnan(occupancy_map_flt(:,:,bbb)))
    xlabel('X (m)')
    ylabel('Y (m)')
    title(sprintf('%s-%s: occupancy (s)',bnames{bb},bnames{bbb}),'Interpreter','none')
    colormap jet
    colorbar

    % Save
    figDir = fullfile(saveDir,'spatial_occupancy');
    if isempty(dir(figDir)); mkdir(figDir); end
    saveas(fig,fullfile(figDir,sprintf('%s.png',bnames{bbb})))
    close(fig)
end