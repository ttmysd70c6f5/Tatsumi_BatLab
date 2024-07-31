function plotVelocityDistribution(SocialFlights,observe_bat,saveDir,savePrefix)

% Behavioral variables
n_bat = length(SocialFlights);

for flying_bat = [1:observe_bat-1,observe_bat+1:n_bat]
    bname = SocialFlights(flying_bat).bname;
    clu_list = SocialFlights(flying_bat).clu_list;
    n_clu = length(clu_list);
    
    %%% FIGURE
    fig = figure('Visible','off');
    % fig = figure;
    set(fig, 'units','normalized','Position',[0 0.2 1 0.6]);
    tl = tiledlayout(2,5,"TileSpacing",'compact');
    title(tl,sprintf('Bat %s',bname))
    for clu = 1:n_clu
        idx_clu_now = SocialFlights(flying_bat).clu_idx == clu_list(clu);
        idx_social = SocialFlights(flying_bat).social;
        idx_asocial = ~SocialFlights(flying_bat).social;
        
        v_social = SocialFlights(flying_bat).v_abs_obs(idx_clu_now & idx_social); % absolute acceleration of observor bat
        v_asocial = SocialFlights(flying_bat).v_abs_obs(idx_clu_now & idx_asocial); % absolute acceleration of observor bat
        n_flight_social = length(v_social);
        n_flight_asocial = length(v_asocial);
        
        nexttile
        hold on
        % social
        histogram(v_social,0:0.2:1,'FaceColor','r','FaceAlpha',0.3)
        histogram(v_asocial,0:0.2:1,'FaceColor','k','FaceAlpha',0.3)
        hold off

        title(sprintf('Clu%d: Social:%d,Non-social:%d',clu_list(clu),n_flight_social,n_flight_asocial));
        xlabel('')
        set(gca,'FontSize',14)
    end
    
    
    if isempty(dir(saveDir)); mkdir(saveDir); end
    saveas(fig,fullfile(saveDir,sprintf('%s%s.png',savePrefix,bname)))
end