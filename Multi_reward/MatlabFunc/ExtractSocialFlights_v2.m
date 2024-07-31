function [SocialFlights, AsocialFlights] = ExtractSocialFlights_v2(BhvData,Flights,TrjCluster,observe_bat,target)
% Extract social/non-social flights from the Flights structure

% Behavioral variables
n_bats = BhvData.n_tags;
bflying = BhvData.bflying; % 1 if at resting sites
blanding = 1 - bflying;
r = BhvData.r;
t = BhvData.t;

% Params
th_social = 0.5; % threshold for social flight (m)
th_fd = 0.3;

% Extraction of social/nonsocial flights
n_flight = Flights(observe_bat).n_flight; % # of flights
t_land_idx = Flights(observe_bat).t_land_idx; % time index of landings

for flying_bat = [1:observe_bat-1,observe_bat+1:n_bats]
    % Identify when the other bats take-off/land while the observing bat is stationary
    t_flight_idx = []; % conspecific land timings
    for ll = 1:n_flight
        t_start = t_land_idx(ll,1); t_end = t_land_idx(ll,2);
        if strcmp(target,'takeoff')
            is_oth_flight = diff([bflying(t_start,flying_bat); bflying(t_start:t_end,flying_bat)]) == 1;
        elseif strcmp(target,'landing')
            is_oth_flight = diff([blanding(t_start,flying_bat); blanding(t_start:t_end,flying_bat)]) == 1;
        end
        t_flight_idx = [t_flight_idx; t_land_idx(ll,1) + find(is_oth_flight) - 1]; % Index of others' takeoffs        
    end
    is_social = vecnorm(diff(r(t_flight_idx,:,[flying_bat,observe_bat]),1,3),2,2) < th_social;
    idx_social = is_social;
    idx_asocial = ~is_social;
    
    % Exclude flights that observe bat is moving
    is_stationary = bflying(t_flight_idx,observe_bat) == 0;
    idx_social = is_social & is_stationary;
    idx_asocial = ~is_social & is_stationary;
        
    % Exclude flights that observe bat is at feeders
    r_fd = [2.77,0.82,1.75; 2.79,-0.99,1.64; 2.78,1.29,0.84; 2.78,-1.43,0.80]; % Feeder position
    [r_fd,~] = r_feeder_calibrate('Y:\users\Tatsumi\Data\Flight_room_test\240120_feeder_calibration');
    is_feed = false(length(t_flight_idx),4);
    for fd = 1:4
        is_feed(:,fd) = vecnorm(r(t_flight_idx,:,observe_bat) - r_fd(fd,:),2,2) < th_fd;
    end
    is_feed = any(is_feed,2); % 1 if on feeders
    idx_social = is_social & is_stationary & ~is_feed;
    idx_asocial = ~is_social & is_stationary & ~is_feed;
        
    SocialFlights(flying_bat).t_idx = t_flight_idx(idx_social);
    SocialFlights(flying_bat).t_sec = t(SocialFlights(flying_bat).t_idx);
    SocialFlights(flying_bat).a_abs_observe = BhvData.ac_abs(SocialFlights(flying_bat).t_idx,observe_bat);
    % SocialFlights(flying_bat).clu_idx = TrjCluster(flying_bat).clu_idx(idx_social);

    AsocialFlights(flying_bat).t_idx = t_flight_idx(idx_asocial);
    AsocialFlights(flying_bat).t_sec = t(AsocialFlights(flying_bat).t_idx);
    AsocialFlights(flying_bat).a_abs_observe = BhvData.ac_abs(AsocialFlights(flying_bat).t_idx,observe_bat);
    % AsocialFlights(flying_bat).clu_idx = TrjCluster(flying_bat).clu_idx(idx_asocial);
end
