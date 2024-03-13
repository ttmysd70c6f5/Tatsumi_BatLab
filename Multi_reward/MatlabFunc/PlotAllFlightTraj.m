function PlotAllFlightTraj(BhvData,RewardData,figDir,varargin)
% BhvData: Extracted behavior data
% figDir: path of figure output directory


p = inputParser;
addRequired(p, 'BhvData');
addRequired(p, 'RewardData');
addRequired(p, 'figDir');
addOptional(p, 'include', 'all'); % all = all landings, fd = only feeder landings

parse(p, BhvData, RewardData, figDir, varargin{:});
extractMode = p.Results.include;

% output directory
if ~isempty(dir(figDir))
    disp('Figure directory already exists.')
    return
else
    mkdir(figDir)
    
    % Detect landings
    switch extractMode
        case 'all'
            Landings = extractLandings(BhvData,RewardData,'include','all');
        case 'feeder'
            Landings = extractLandings(BhvData,RewardData,'include','feeder');
    end

    t = BhvData.t; % timestamps
    T = BhvData.T; % Length of recording in CDP
    r_lim = BhvData.r_lim; % dimension of room
    blanding = 1 - BhvData.bflying; % 1 if landing
    n_tags = BhvData.n_tags; % number of animals
    r = BhvData.r; % position of bats
    
    % Feeder coordinates
    % r_fd = [2.77,0.82,1.75; 2.79,-0.99,1.64; 2.78,1.29,0.84; 2.78,-1.43,0.80]; % Feeder position
    [r_fd,~] = r_feeder_calibrate('Y:\users\Tatsumi\Data\Flight_room_test\240120_feeder_calibration');
    
    for bb = 1:n_tags
    
        t_land = Landings(bb).t_land_sec;
        n_land = Landings(bb).n_land;
        is_fd_land = Landings(bb).is_fd_land;
    
        range_draw = 20; % range to visualize around the feeder activation timing
        pre_t = 5; % time to draw before landing starts
        
        figure('Units','normalized','Position',[0.15,0.05,0.7,0.9],'Visible','off')
        tl = tiledlayout(5,5,"TileSpacing",'loose');
        title(tl,sprintf('Bat%d: %d sec after landings',bb,range_draw))
        fig_cnt = 1;
        save_cnt = 1;
        for ll = 1:n_land
            idx_t = find(ismember(round(t,2),round(t_land(ll),2)));
                if idx_t > 100*pre_t && idx_t < T - 100*range_draw
                    r_plot = r(idx_t-100*pre_t:idx_t+100*range_draw,:,bb);
                    t_plot = (-1*pre_t:0.01:range_draw)' + t(idx_t);
                    l_plot = blanding(idx_t-100*pre_t:idx_t+100*range_draw,bb);
                elseif idx_t <= 100*pre_t
                    r_plot = r(1:idx_t+100*range_draw,:,bb);
                    t_plot = (-1*(idx_t-1)/100:0.01:range_draw)' + t(idx_t);
                    l_plot = blanding(1:idx_t+100*range_draw,bb);
                else
                    r_plot = r(idx_t-100*pre_t:end,:,bb);
                    t_plot = (-1*pre_t:0.01:(T-idx_t)/100)' + t(idx_t);
                    l_plot = blanding(idx_t-100*pre_t:end,bb);
                end
            n = length(t_plot);
            CD = [uint8(jet(n)*255) uint8(ones(n,1))].';
        
            % t_plot = t(idx_feeder-100*60:idx_feeder+100*30);
        
            nexttile([1,4])
            hold on
            plot(t_plot,r_plot,'LineWidth',1);
            plot([t(idx_t),t(idx_t)],[-3,3],'k')
            area(t_plot,l_plot,'FaceAlpha',0.2,'LineStyle','none','FaceColor',[0.3010 0.7450 0.9330])
            hold off
            xlim([t(idx_t)-pre_t,t(idx_t)+range_draw])
            ylim([-3,3])
            if is_fd_land(ll)
                title(sprintf('%dth landing: feeder',ll))
            else
                title(sprintf('%dth landing: non-feeder',ll))
            end
            xlabel('Time (sec)')
            ylabel('Position (m)')
            legend({'x','y','z','Beam blocked'})
        
            nexttile([1,1])
            hold on
            p = plot(r_plot(:,1),r_plot(:,2));
            drawnow
            set(p.Edge, 'ColorBinding','interpolated', 'ColorData',CD)
            scatter(r(idx_t,1,bb),r(idx_t,2,bb),20,'k','filled')
    
            rectangle('Position',[r_fd(1,1)-0.1,r_fd(1,2)-0.1,0.2,0.2],'EdgeColor','k')
            rectangle('Position',[r_fd(2,1)-0.1,r_fd(2,2)-0.1,0.2,0.2],'EdgeColor','k')
            rectangle('Position',[r_fd(3,1)-0.1,r_fd(3,2)-0.1,0.2,0.2],'EdgeColor','c')
            rectangle('Position',[r_fd(4,1)-0.1,r_fd(4,2)-0.1,0.2,0.2],'EdgeColor','c')
            hold off
            title('Top view')
            xlabel('x (m)')
            ylabel('y (m)')
            xlim([r_lim(1,1),r_lim(1,2)])
            ylim([r_lim(2,1),r_lim(2,2)])
        
            % nexttile([1,1])
            % hold on
            % p = plot(r_plot(:,2),r_plot(:,3));
            % drawnow
            % set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
            % scatter(r(idx_feeder,2,bb),r(idx_feeder,3,bb),20,'k','filled')
            % hold off
            % xlabel('y')
            % ylabel('z')
            % xlim([r_lim(2,1),r_lim(2,2)])
            % ylim([r_lim(3,1),r_lim(3,2)])
        
            if fig_cnt < 5
                fig_cnt = fig_cnt + 1;
            else
                if ll ~= n_land
                    saveas(gcf,fullfile(figDir,['AllLandings_bat' num2str(bb) '_' num2str(save_cnt) '.png']));
                    save_cnt = save_cnt + 1;
                    close all
    
                    figure('Units','normalized','Position',[0.15,0.05,0.7,0.9],'Visible','off')
                    tl = tiledlayout(5,5,"TileSpacing",'loose');
                    title(tl,sprintf('Bat%s',num2str(bb)))
                    fig_cnt = 1;
                end
            end
        end
        
        saveas(gcf,fullfile(figDir,['AllLandings_bat' num2str(bb) '_' num2str(save_cnt) '.png']));
        close all
    end
end

end