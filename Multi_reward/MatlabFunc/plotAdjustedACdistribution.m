function plotAdjustedACdistribution(SocialFlights,observe_bat,saveDir,savePrefix)

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
        % idx_clu_now = SocialFlights(flying_bat).clu_idx == clu_list(clu);
        % idx_social = SocialFlights(flying_bat).social & ~SocialFlights(flying_bat).obs_at_fd & SocialFlights(flying_bat).stationary & SocialFlights(flying_bat).is_include;
        % idx_asocial = ~SocialFlights(flying_bat).social & ~SocialFlights(flying_bat).obs_at_fd & SocialFlights(flying_bat).stationary & SocialFlights(flying_bat).is_include;
        % 
        ac_social = SocialFlights(flying_bat).a_adjusted{clu,1}; % absolute acceleration of observor bat
        ac_asocial = SocialFlights(flying_bat).a_adjusted{clu,2}; % absolute acceleration of observor bat
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