function RasterPlot_clusteredFlights(sp,ephys_idx,bhv_idx,figDir)
bhv_idx
bname_bhv = sp(bhv_idx).batname; % bat name for behavior data
bname_ephys = sp(ephys_idx).batname; % bat name for ephys data

figDir = fullfile(figDir,bname_bhv);
if isempty(dir(figDir))
    mkdir(figDir)
end

n_unit = length(sp(bhv_idx).sptimes(1).sptimes_flight);
n_clu = length(unique(sp(bhv_idx).trj_cluster));
clu_list = unique(sp(bhv_idx).trj_cluster);

for unit = 1:n_unit
    % unit = 1;
    unit_depth = sp(1).sptimes(1).unit_depth(unit).unit_depth; % depth of SU/MUA
            
    % nexttile; hold on
    f1 = figure('Visible','off'); hold on
    set(f1, 'units','normalized','Position',[0.2 0.2 0.6 0.6]);
    tl = tiledlayout(6,9,"TileSpacing",'compact');
    % title(sprintf('Unit %d',unit))
    title(tl,sprintf('Flight #%s , ephys #%s: SU %d, depth %dum',bname_bhv,bname_ephys,unit,unit_depth))
    
    for clu = 1:n_clu
        clu_name = clu_list(clu);

        clu_idx = sp(bhv_idx).trj_cluster == clu_name; % extract the flight cluster
        
        clu_size = sum(clu_idx);
        sessions = sp(bhv_idx).session(clu_idx);
        len_flight = sp(bhv_idx).len_flight_sec(clu_idx);
        len_land = sp(bhv_idx).len_land_sec(clu_idx);
        sptimes_clu = sp(bhv_idx).sptimes(clu_idx);
        r_flight = sp(bhv_idx).r_flight(clu_idx);

        for session = 1:length(unique(sessions))
            sptimes_clu_session = sptimes_clu(sessions==session);
            [sorted_flight_session,sorted_idx] = sort(len_flight(sessions==session),'ascend');
            len_flight(sessions==session) = sorted_flight_session;
            sptimes_clu(sessions==session) = sptimes_clu_session(sorted_idx);
        end
    
        % figure('Visible','off')
        % set(gcf, 'units','normalized','Position',[0 0 1 1]);
        % tl = tiledlayout(5,10,"TileSpacing",'compact');
        % title(tl,sprintf('Raster plot during flight: bat 14662, cluster %d',clu_name))
        
        

        % Plot trajectory
        nexttile; hold on
        for trj = 1:clu_size
            x = r_flight{trj}(:,1,bhv_idx);
            y = r_flight{trj}(:,2,bhv_idx);
            z = r_flight{trj}(:,3,bhv_idx);
            n = length(x);
            GradientColor = [uint8(jet(n)*255) uint8(ones(n,1))].';
            p = plot3(x,y,z);
            drawnow
            set(p.Edge, 'ColorBinding','interpolated', 'ColorData',GradientColor)
        end
        hold off
        title(sprintf('Cluster %d',clu_name))
        xlabel('x'); ylabel('y'); zlabel('z')
        xlim([-3,3]); ylim([-3,3]); zlim([0,2.5])
        view([0,90])

        % Plot raster plots
        nexttile
        for i = 1:clu_size
            sptimes = sptimes_clu(i).sptimes_flight(unit).sptimes_usec / 1e6;
            for ii = 1:length(sptimes)
                line([sptimes(ii) sptimes(ii)], [i-1 i],'Color','k')
            end
            line([len_flight(i) len_flight(i)],[i-1 i],'Color','r','LineWidth',1)
        end

        session_list = unique(sessions);
        for ss = 1:length(session_list)
            session = session_list(ss);
            if mod(session,2) == 0
                session_start_idx = find(sessions==session,1,'first');
                session_end_idx = find(sessions==session,1,'last');
                yregion(session_start_idx-1,session_end_idx,'FaceColor',[80,190,225]/255)
            end
        end
        hold off
        xlabel('time (sec)')
        title(sprintf('Cluster %d',clu_name))
        % ylabel(bhv_type)

        % PETH
        nexttile
        edges = [0:0.1:5]; % edges
        peth = zeros(1,length(edges)-1); % initialize the PETH
        for i = 1:clu_size
            sptimes = sptimes_clu(i).sptimes_flight(unit).sptimes_usec / 1e6; % spike times in sec
            [N,~] = histcounts(sptimes,edges); % count the spike times in each bin
            peth = peth + N; % Add current flight
        end
        % peth(i,:) = normalize(peth(i,:),'zscore'); % normalize peth
        peth = peth ./ clu_size;
        bar([0.05:0.1:4.95],peth) % Plot PETH as a bar graph
        xlim([0,5])
        xlabel('Time (sec)')
        ylabel('z-scored firing rate')
        
    end

    saveas(f1,fullfile(figDir,sprintf('raster_flight_%s_ephys_%s_SU%d.jpg',bname_bhv,bname_ephys,unit)))
    close(f1)
    % saveas(gcf,fullfile(figDir,sprintf('raster_ephys_14662_flight_bat%s_cluster%d.jpg',bname,clu_name)))
    % close gcf
end

end