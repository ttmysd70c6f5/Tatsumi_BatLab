function plotACdistribution(str,obs,saveDir,savePrefix)

% Behavioral variables
n_bat = length(str);

for flying_bat = [1:obs-1,obs+1:n_bat]
    bname = str(flying_bat).bname;
    clu_list = str(flying_bat).clu_list;
    n_clu = length(clu_list);
    
    %%% FIGURE
    fig = figure('Visible','off');
    % fig = figure;
    set(fig, 'units','normalized','Position',[0 0.2 1 0.6]);
    tl = tiledlayout(2,5,"TileSpacing",'compact');
    title(tl,sprintf('Bat %s',bname))

    for clu = 1:n_clu
        idx_clu_now = str(flying_bat).clu_idx == clu_list(clu);
        idx_social = str(flying_bat).social;
        idx_asocial = ~str(flying_bat).social;
        
        ac_social = str(flying_bat).a_abs_obs(idx_clu_now & idx_social); % absolute acceleration of observor bat
        ac_asocial = str(flying_bat).a_abs_obs(idx_clu_now & idx_asocial); % absolute acceleration of observor bat
        n_flight_social = length(ac_social);
        n_flight_asocial = length(ac_asocial);
        
        nexttile
        hold on
        % social
        histogram(ac_social,0:0.25:3,'FaceColor','r','FaceAlpha',0.3)
        histogram(ac_asocial,0:0.25:3,'FaceColor','k','FaceAlpha',0.3)
        hold off

        title(sprintf('Clu%d: Social:%d,Non-social:%d',clu_list(clu),n_flight_social,n_flight_asocial));
        xlabel('')
        set(gca,'FontSize',14)
    end
    
    
    if isempty(dir(saveDir)); mkdir(saveDir); end
    saveas(fig,fullfile(saveDir,sprintf('%s%s.png',savePrefix,bname)))
end