function RasterPlot_unclusteredBehavior(sp,ephys_bat,bat_idx,bhv_types)
% Generate the raster plot during the landings/flights of the bat specified by bat_idx.
% INPUT
% sp: the struct containing the spike times during the flight/flights
% ephys_bat: the ID of the ephys bat
% bat_idx: the ID of the bat to be extracted for the behavior
% bhv_types: if bhv_types is 'landing', the raster plot is generated for
% the spike times during each landing. On the other hand, the raster plot
% is generated for the spike times during each flight if 'flight' is given
% as bhv_types.

for bhv = 1:length(bhv_types)
    bhv_type = bhv_types{bhv};

    % Make figure directory
    figDir = fullfile(pwd,'analysis','figure',sprintf('RasterPlot_%s_Ephys%s_nonclustering',bhv_type,ephys_bat));
    if isempty(dir(figDir))
        mkdir(figDir)
    end

    % Run for all bats
    for b = 1:length(bat_idx)
        bb = bat_idx(b); % index of bat
        bname = sp(bb).batname; % bat name
        n_unit = length(sp(bb).sptimes(1).sptimes_flight); % number of single units
        
        % Run for all single units
        for unit = 1:n_unit
            unit_depth = sp(1).sptimes(1).unit_depth(unit).unit_depth; % depth of SU/MUA

            % Figure settings
            f1 = figure('Visible','off');
            set(f1, 'units','normalized','Position',[0.2 0.2 0.8 0.6]);
            tl = tiledlayout(1,3,"TileSpacing",'compact');
            title(tl,sprintf('%s: ephys %s, SU %d, depth %dum',bhv_type, ephys_bat, unit, unit_depth))

            % Run for all landing spots
            for feeder_id = 0:2
                switch bhv_type
                    case 'landing'
                        len_bhv_sec = sp(bb).len_land_sec; % length of landing in sec
                    case 'flight'
                        len_bhv_sec = sp(bb).len_flight_sec; % landing of flight in sec
                end
    
                % Extract the landings/flights on the selected feeder/non-feeder
                if feeder_id ~= 0
                    bhv_idx = sp(bb).is_fd_land & sp(bb).closest_feeder(:,2) == feeder_id & len_bhv_sec < 90; % logical index of landings for the feeder
                else
                    bhv_idx = ~sp(bb).is_fd_land & len_bhv_sec < 90; % landings on non-feeder
                end
            
                num_bhv = sum(bhv_idx); % number of flights/landings
                sessions = sp(bb).session(bhv_idx); % current session
                sptimes_bhv = sp(bb).sptimes(bhv_idx); % spike times during flights/landings
                len_bhv_sec =  len_bhv_sec(bhv_idx); % length of flights/landings in sec
            
                % Sort the data according to the length of landing
                % Sort for each session
                for session = 1:length(unique(sessions))
                    sptimes_bhv_session = sptimes_bhv(sessions==session); % spike times for the current session
                    [sorted_bhv_session,sorted_idx] = sort(len_bhv_sec(sessions==session),'ascend'); % sort based on the length of flights/landings
                    len_bhv_sec(sessions==session) = sorted_bhv_session; % Replace the old values
                    sptimes_bhv(sessions==session) = sptimes_bhv_session(sorted_idx); % replace the old values
                end

                % Move to the new axis
                nexttile; hold on
                title(sprintf('Feeder %d',feeder_id))

                % Run for all trials
                for i = 1:num_bhv
                    switch bhv_type
                        case 'landing'
                            sptimes = sptimes_bhv(i).sptimes_land(unit).sptimes_usec / 1e6; % spike times in usec
                        case 'flight'
                            sptimes = sptimes_bhv(i).sptimes_flight(unit).sptimes_usec / 1e6; % spike times in usec
                    end
                    
                    % Draw the raster plot
                    % switch bhv_type
                    %     case 'landing'
                    %         for ii = 1:length(sptimes)
                    %             line([-sptimes(ii) -sptimes(ii)], [i-1 i],'Color','k')
                    %         end
                    %         line([-len_bhv_sec(i) -len_bhv_sec(i)],[i-1 i],'Color','r','LineWidth',1)
                    %         drawnow
                    %     case 'flight'
                    %         for ii = 1:length(sptimes)
                    %             line([sptimes(ii) sptimes(ii)], [i-1 i],'Color','k')
                    %         end
                    %         line([len_bhv_sec(i) len_bhv_sec(i)],[i-1 i],'Color','r','LineWidth',1)
                    %         drawnow
                    % end


                    for ii = 1:length(sptimes)
                        line([sptimes(ii) sptimes(ii)], [i-1 i],'Color','k')
                    end
                    line([len_bhv_sec(i) len_bhv_sec(i)],[i-1 i],'Color','r','LineWidth',1)
                    drawnow

                    
                end
        
                % Change the background color for even sessions
                for session = 1:length(unique(sessions))
                    if mod(session,2) == 0
                        session_start_idx = find(sessions==session,1,'first');
                        session_end_idx = find(sessions==session,1,'last');
                        yregion(session_start_idx-1,session_end_idx,'FaceColor',[80,190,225]/255)
                    end
                end
                hold off
            end

            % Save figure
           saveas(f1,fullfile(figDir,sprintf('raster_%s_%s_ephys_%s_SU%d.jpg',bname,bhv_type,ephys_bat,unit)))
           close(f1)
        end
    end
end

end