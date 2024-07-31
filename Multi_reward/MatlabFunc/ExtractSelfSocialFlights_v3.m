function out = ExtractSelfSocialFlights_v3(BhvData,Flights,TrjCluster,fly_bat,target,t_start,t_end)
% Extract social/non-social flights from the Flights structure
% Does MATTER if the observing bat is stationary or flying

% Behavioral variables
n_bat = BhvData.n_tags;
bflying = BhvData.bflying;
r = BhvData.r;
t = BhvData.t;

% Params
th_social = 0.5; % threshold for social flight (m)
th_fd = 0.3;

% Extraction
% observe_bat = find(strcmp(mtData.bname_ephys,mtData.bname)); % Ephys bat

for target_bat = [1:fly_bat-1,fly_bat+1:n_bat] % Flying bat
    % Identify social flights
    if strcmp(target,'takeoff')
        t_flight_idx = Flights(fly_bat).t_flight_idx(:,1); % takeoff time index
        r_extracted = cell2mat(cellfun(@(x) x(1,:,[fly_bat,target_bat]), Flights(fly_bat).r_flight, 'UniformOutput',false)); % position of the flying/observing bat at takeoffs by the flying bat
    elseif strcmp(target, 'landing')
        t_flight_idx = Flights(fly_bat).t_land_idx(:,1); % landing time index
        r_extracted = cell2mat(cellfun(@(x) x(1,:,[fly_bat,target_bat]), Flights(fly_bat).r_land, 'UniformOutput',false)); % position of the flying/observing bat at takeoffs by the flying bat
    end
        
    d = vecnorm(diff(r_extracted,1,3),2,2); % distance between flying/observing bats at takeoffs
    is_social = d < th_social;

    % Exclude flights that observe bat is moving
    is_stationary = bflying(t_flight_idx,target_bat) == 0;

    % Exclude flights that observe bat is at feeders
    % r_fd = [2.77,0.82,1.75; 2.79,-0.99,1.64; 2.78,1.29,0.84; 2.78,-1.43,0.80]; % Feeder position
    [r_fd,~] = r_feeder_calibrate('Y:\users\Tatsumi\Data\Flight_room_test\240120_feeder_calibration');
    is_feed = false(length(t_flight_idx),4);
    for fd = 1:4
        is_feed(:,fd) = vecnorm(r(t_flight_idx,:,target_bat) - r_fd(fd,:),2,2) < th_fd;
    end
    is_feed = any(is_feed,2); % 1 if on feeders

    % Exclude the flight at very biggining or end
    t_sec = t(t_flight_idx);
    is_include = t_sec > t_start & t_sec < t_end;

    idx_to_save = is_stationary & ~is_feed & is_include;

    % Save to output structure
    out(target_bat).bname = Flights(target_bat).bname;
    out(target_bat).target = target;
    out(target_bat).social = is_social(idx_to_save);
    % SocialFlights(flying_bat).stationary = is_stationary;
    % SocialFlights(flying_bat).obs_at_fd = is_feed;
    
    if strcmp(target,'takeoff')
        out(target_bat).t_idx = Flights(fly_bat).t_flight_idx(idx_to_save,1);
    elseif strcmp(target, 'landing')
        out(target_bat).t_idx = Flights(fly_bat).t_land_idx(idx_to_save,2);
    end
    out(target_bat).t_sec = t(out(target_bat).t_idx);
    out(target_bat).a_abs_obs = BhvData.ac_abs(out(target_bat).t_idx,target_bat);
    out(target_bat).v_abs_obs = BhvData.v_abs(out(target_bat).t_idx,target_bat);
    out(target_bat).r_obs = r(out(target_bat).t_idx,:,target_bat);
    out(target_bat).r_fly = r(out(target_bat).t_idx,:,fly_bat);
    out(target_bat).d = d(idx_to_save);
    out(target_bat).clu_idx = TrjCluster(fly_bat).clu_idx(idx_to_save);
    out(target_bat).clu_list = TrjCluster(fly_bat).clu_list;
    out(target_bat).n_flight = length(out(target_bat).social);
end