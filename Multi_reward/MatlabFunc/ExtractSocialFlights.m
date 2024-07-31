function [SocialFlights, AsocialFlights] = ExtractSocialFlights(BhvData,Flights,observe_bat,target)
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
SocialFlights = struct([]);
AsocialFlights = struct([]);
for flying_bat = [1:observe_bat-1,observe_bat+1:n_bat] % Flying bat
    % Identify social flights
    if strcmp(target,'takeoff')
        t_flight_idx = Flights(flying_bat).t_flight_idx(:,1); % takeoff time index
        r_extracted = cell2mat(cellfun(@(x) x(1,:,[observe_bat,flying_bat]), Flights(flying_bat).r_flight, 'UniformOutput',false)); % position of the flying/observing bat at takeoffs by the flying bat
    elseif strcmp(target, 'landing')
        t_flight_idx = Flights(flying_bat).t_land_idx(:,1); % landing time index
        r_extracted = cell2mat(cellfun(@(x) x(1,:,[observe_bat,flying_bat]), Flights(flying_bat).r_land, 'UniformOutput',false)); % position of the flying/observing bat at takeoffs by the flying bat
    end
        
    d = vecnorm(diff(r_extracted,1,3),2,2); % distance between flying/observing bats at takeoffs
    is_social = d < th_social;

    % Exclude flights that observe bat is moving
    is_stationary = bflying(t_flight_idx,observe_bat) == 0;

    % Exclude flights that observe bat is at feeders
    % r_fd = [2.77,0.82,1.75; 2.79,-0.99,1.64; 2.78,1.29,0.84; 2.78,-1.43,0.80]; % Feeder position
    [r_fd,~] = r_feeder_calibrate('Y:\users\Tatsumi\Data\Flight_room_test\240120_feeder_calibration');
    is_feed = false(length(t_flight_idx),4);
    for fd = 1:4
        is_feed(:,fd) = vecnorm(r(t_flight_idx,:,observe_bat) - r_fd(fd,:),2,2) < th_fd;
    end
    is_feed = any(is_feed,2); % 1 if on feeders

    if strcmp(target,'takeoff')
        SocialFlights(flying_bat).t_idx = Flights(flying_bat).t_flight_idx(is_social & is_stationary & ~is_feed,:);
        AsocialFlights(flying_bat).t_idx = Flights(flying_bat).t_flight_idx(~is_social & is_stationary & ~is_feed,:);
    elseif strcmp(type, 'landing')
        SocialFlights(flying_bat).t_idx = Flights(flying_bat).t_land_idx(is_social & is_stationary & ~is_feed,:);
        AsocialFlights(flying_bat).t_idx = Flights(flying_bat).t_land_idx(~is_social & is_stationary & ~is_feed,:);
    end
    SocialFlights(flying_bat).t_sec = t(SocialFlights(flying_bat).t_idx);
    SocialFlights(flying_bat).a_abs_observe = BhvData.ac_abs(SocialFlights(flying_bat).t_idx,observe_bat);  
    AsocialFlights(flying_bat).t_sec = t(AsocialFlights(flying_bat).t_idx);    
    AsocialFlights(flying_bat).a_abs_observe = BhvData.ac_abs(AsocialFlights(flying_bat).t_idx,observe_bat);
end