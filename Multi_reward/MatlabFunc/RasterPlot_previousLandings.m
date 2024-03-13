function RasterPlot_previousLandings(sp,ephys_bat,bat_idx)

bhv_type = 'landing';

figDir = fullfile(pwd,'analysis','figure',sprintf('RasterPlot_%s_before_flight2feeder_Ephys%s',bhv_type,ephys_bat));
if isempty(dir(figDir))
    mkdir(figDir)
end

for b = 1:length(bat_idx)
    bb = bat_idx(b);
    bname = sp(bb).batname; % bat name

               
    % Generate raster plot for each unit for all flight clusters
    n_unit = length(sp(bb).sptimes(1).sptimes_flight);
    for unit = 1:n_unit
        unit_depth = sp(1).sptimes(1).unit_depth(unit).unit_depth; % depth of SU/MUA

        % Figure settings
        % f1 = figure('Visible','off');
        f1 = figure;
        set(gcf, 'units','normalized','Position',[0.2 0.2 0.8 0.6]);
        tl = tiledlayout(1,3,"TileSpacing",'compact');
        title(tl,sprintf('Raster plot during %s: ephys %s, SU %d, depth %dum',bhv_type, ephys_bat, unit, unit_depth))

        for feeder_id = 0:2
        % for feeder_id = 2
            switch bhv_type
                case 'landing'
                    len_bhv_sec = sp(bb).len_land_sec;
                case 'flight'
                    len_bhv_sec = sp(bb).len_flight_sec;
            end

            % Extract landings on the feeder or non-feeder places
            if feeder_id ~= 0
                bhv_idx = sp(bb).is_fd_land & sp(bb).closest_feeder(:,2) == feeder_id & [0;len_bhv_sec(1:end-1)] < 70; % logical index of landings for the feeder
            else
                bhv_idx = ~sp(bb).is_fd_land & [0;len_bhv_sec(1:end-1)] < 70; % landings on non-feeder
            end
        
            if find(bhv_idx,1,'first')==1
                bhv_idx(find(bhv_idx,1,'first')) = false; % exclude the first landing
            end
            bhv_idx = [bhv_idx(2:end);false]; % extract the previous landing
                
            num_bhv = sum(bhv_idx);
            sessions = sp(bb).session(bhv_idx);
            sptimes_bhv = sp(bb).sptimes(bhv_idx);
            len_bhv_sec =  len_bhv_sec(bhv_idx);
        
            % Sort the data according to the length of landing in each session
            for session = 1:length(unique(sessions))
                sptimes_bhv_session = sptimes_bhv(sessions==session);
                [sorted_bhv_session,sorted_idx] = sort(len_bhv_sec(sessions==session),'ascend');
                len_bhv_sec(sessions==session) = sorted_bhv_session;
                sptimes_bhv(sessions==session) = sptimes_bhv_session(sorted_idx);
            end
            
            nexttile; hold on
            title(sprintf('Feeder %d',feeder_id))
            for i = 1:num_bhv
                switch bhv_type
                    case 'landing'
                        sptimes = sptimes_bhv(i).sptimes_land(unit).sptimes_usec / 1e6;
                    case 'flight'
                        sptimes = sptimes_bhv(i).sptimes_flight(unit).sptimes_usec / 1e6;
                end
                
                for ii = 1:length(sptimes)
                    line([-sptimes(ii) -sptimes(ii)], [i-1 i],'Color','k')
                end
                line([-len_bhv_sec(i) -len_bhv_sec(i)],[i-1 i],'Color','r','LineWidth',1)
                drawnow
            end
    
            for session = 1:length(unique(sessions))
                if mod(session,2) == 0
                    session_start_idx = find(sessions==session,1,'first');
                    session_end_idx = find(sessions==session,1,'last');
                    yregion(session_start_idx-1,session_end_idx,'FaceColor',[80,190,225]/255)
                end
            end
            hold off
            xlim([-80,0])
        end
       saveas(f1,fullfile(figDir,sprintf('raster_ephys_%s_%s_bat%s_unit%d.jpg',bhv_type,ephys_bat,bname,unit)))
       close(f1)
    end
end
end