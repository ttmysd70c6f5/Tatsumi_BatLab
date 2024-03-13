function sp = extractSpikeDuringBehavior(spikeData,Landings,TrjCluster,mtData,unit_label)
% Align the neural and behavioral data into a struct.
% INPUT
% spikeData: spike times extracted by loadSpikeData
% Landings: landing/flight extracted by extractLandings
% TrjCluster: clustered trajectories by TrajClusteringDBSCAN
% mtData: experiment info
% unit_label: su = single unit only, mua = multi-unit only

sp = struct([]); % spike activities during flight/landing
for bb = 1:length(Landings)
    % bb = 1;
    sp(bb).batname = mtData.bname{bb};
    sp(bb).t_land_sec = Landings(bb).t_land_sec;
    sp(bb).t_flight_sec = Landings(bb).t_flight_sec;
    sp(bb).len_land_sec = sp(bb).t_land_sec(:,2) - sp(bb).t_land_sec(:,1); % length of flight/landing
    sp(bb).len_flight_sec = sp(bb).t_flight_sec(:,2) - sp(bb).t_flight_sec(:,1); % length of flight/landing
    sp(bb).r_flight = Landings(bb).r_flight;
    sp(bb).r_land = Landings(bb).r_land;
    sp(bb).is_fd_land = Landings(bb).is_fd_land;
    sp(bb).closest_feeder = Landings(bb).closest_feeder;
    sp(bb).session = Landings(bb).session_land;
    sp(bb).trj_cluster = TrjCluster(bb).cluster;
    
    
    n_land = Landings(bb).n_land;
    
    for ll = 1:n_land
        % ll = 1;
        t_land_sec = sp(bb).t_land_sec(ll,:);
        t_flight_sec = sp(bb).t_flight_sec(ll,:);
        
        switch unit_label
            case 'su'
                n_unit = spikeData.numGoodUnits;
                % sptimes = struct([]); % spike times
                for unit = 1:n_unit
                    % unit = 1;
                    unit_depth = spikeData.good_units.depth(unit); % depth of the unit
                    sptimes_unit = spikeData.good_units.spikeTimes_usec{unit}; % spike times of the current unit
                    sptimes_unit_land = sptimes_unit(sptimes_unit >= t_land_sec(1)*1e6 & sptimes_unit <= t_land_sec(2)*1e6) - t_land_sec(1)*1e6; % logical index for the spike times for the current flight/landing in the current unit
                    sptimes_unit_flight = sptimes_unit(sptimes_unit >= t_flight_sec(1)*1e6 & sptimes_unit <= t_flight_sec(2)*1e6) - t_flight_sec(1)*1e6; % logical index for the spike times for the current flight/landing in the current unit
                    
                    sp(bb).sptimes(ll).unit_depth(unit).unit_depth = unit_depth;
                    sp(bb).sptimes(ll).sptimes_land(unit).sptimes_usec = sptimes_unit_land;
                    sp(bb).sptimes(ll).sptimes_flight(unit).sptimes_usec = sptimes_unit_flight;
                end
            case 'mua'
                n_unit = spikeData.numMuaUnits;
                % sptimes = struct([]); % spike times
                for unit = 1:n_unit
                    % unit = 1;
                    unit_depth = spikeData.mua_units.depth(unit); % depth of the unit
                    sptimes_unit = spikeData.mua_units.spikeTimes_usec{unit}; % spike times of the current unit
                    sptimes_unit_land = sptimes_unit(sptimes_unit >= t_land_sec(1)*1e6 & sptimes_unit <= t_land_sec(2)*1e6) - t_land_sec(1)*1e6; % logical index for the spike times for the current flight/landing in the current unit
                    sptimes_unit_flight = sptimes_unit(sptimes_unit >= t_flight_sec(1)*1e6 & sptimes_unit <= t_flight_sec(2)*1e6) - t_flight_sec(1)*1e6; % logical index for the spike times for the current flight/landing in the current unit
                    
                    sp(bb).sptimes(ll).unit_depth(unit).unit_depth = unit_depth;
                    sp(bb).sptimes(ll).sptimes_land(unit).sptimes_usec = sptimes_unit_land;
                    sp(bb).sptimes(ll).sptimes_flight(unit).sptimes_usec = sptimes_unit_flight;
                end
        end

    end
end

end