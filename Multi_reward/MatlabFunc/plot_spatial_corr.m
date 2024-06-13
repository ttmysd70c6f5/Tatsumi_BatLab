function plot_spatial_corr(occupancy_map,spike_map,spikeData,Layers,trjName,bname_ephys)
% Smooth spike and occupancy map
flt_sigma = 1.5;
occupancy_th = 0.2; % threshold for spatial occupancy
spike_map_flt = imgaussfilt(spike_map,flt_sigma);
occupancy_map_flt = imgaussfilt(occupancy_map,flt_sigma);
occupancy_map_flt(occupancy_map < occupancy_th) = nan;
% Compute firing rate map
rate_map_flt = spike_map_flt ./ occupancy_map_flt;


bin_size = 0.15;

n_unit = spikeData.numGoodUnits;
unit_ids = spikeData.good_units.cluster_id;

% unit = 6;
for unit = 1:n_unit

    unit_id = unit_ids(unit);
    depth = spikeData.good_units.depth(unit); % depth of unit
    layers = Layers(strcmp({Layers.trjName},trjName) & strcmp({Layers.batid},bname_ephys)); % specify the layer
    region = 'unknown';
    for i = 1:length(layers.layer)
        if depth > layers.layer(i).depth(1) && depth < layers.layer(i).depth(2)
            region = layers.layer(i).name; % brain region of the current unit
        end
    end
    
    m = size(rate_map_flt,1);
    n = size(rate_map_flt,2);
    SpatialCorr = SpatialAutocorr(rate_map_flt(:,:,unit)); % compute spatial autocorrelation

    
    % fig = figure;
    fig = figure('Visible','off');
    fontsize(fig, 20, "points")
    set(fig, 'units','normalized','Position',[0.2 0.2 0.3 0.4]);
    row_edges = [-1*bin_size*(m-1), bin_size*(m-1)];
    col_edges = [-1*bin_size*(n-1), bin_size*(n-1)];
    h = imagesc(col_edges,row_edges,SpatialCorr);
    set(gca,'YDir','normal','FontSize',16)
    set(h,'AlphaData',~isnan(SpatialCorr))
    xlabel('X (m)')
    ylabel('Y (m)')
    title(sprintf('%s, %s, unit %d',bname_ephys,region,unit_id),'Interpreter','none')
    colormap jet
    colorbar

    saveDir = fullfile(pwd,'analysis','figure','spatial_autocorr',bname_ephys,trjName);
    if isempty(dir(saveDir)); mkdir(saveDir); end
    saveas(fig,fullfile(saveDir,sprintf('spatial_autocorr_Unit%d.jpg',unit_id)))
    close(fig)
end