function social_vector_coding(BhvData,Flights,spikeData,mtData,saveDir)
%% Object coding
% Ephys bat
bb = find(strcmp(mtData.bname_ephys,mtData.bname));

% Load neural and behavioral features
n_bat = BhvData.n_tags;
t = BhvData.t;
r = BhvData.r;
n_unit = spikeData.numGoodUnits;


% Binning settings
bin_size = 0.3; % bin edge length (m)
Xedges = -6:bin_size:6;
Yedges = -6:bin_size:6;
flt_sigma = 1.5; % Sigma for gaussian filter
occupancy_th = 0.1; % threshold for spatial occupancy

% Initialize output structures
occupancy_map_flt = zeros(length(Yedges)-1,length(Xedges)-1,n_bat);
rate_map = zeros(length(Yedges)-1,length(Xedges)-1,n_unit,n_bat);

for bbb = [1:bb-1,bb+1:n_bat]
    fprintf('Bat %d/%d...\n',bbb,n_bat)
    land_idx = true(Flights(bbb).n_flight,1); land_idx(1:min([5,Flights(bbb).n_flight])) = false; % idx of landings considered
    is_fd_land = Flights(bbb).is_fd_land;
    r_land = Flights(bbb).r_land; r_land = r_land(land_idx&~is_fd_land);
    t_land_sec = Flights(bbb).t_land_sec(land_idx&~is_fd_land,:);
    n_land = length(r_land);
    
    % Compute occupancy map of ephys bat while target bat is not moving
    occupancy_map = zeros(length(Yedges)-1,length(Xedges)-1);
    for ll = 1:n_land
        X = r_land{ll}(:,1,bb)-r_land{ll}(:,1,bbb); % Allocentric position
        Y = r_land{ll}(:,2,bb)-r_land{ll}(:,2,bbb);
        [N,~,~] = histcounts2(Y,X,Yedges,Xedges); % row = Y, col = X
        occupancy_map = occupancy_map + N; % append to the occupancy map
    end
    occupancy_map = occupancy_map * 0.1; % convert to sec
    occupancy_map_temp = imgaussfilt(occupancy_map,flt_sigma);
    occupancy_map_temp(occupancy_map < occupancy_th) = nan;
    occupancy_map_flt(:,:,bbb) = occupancy_map_temp;
    
    % Compute firing at allocentric coordinates  
    for unit = 1:n_unit
        if mod(unit,10)==1
            fprintf('Unit %d/%d...\n',unit,n_unit)
        end
        % Initialize output matrices
        spike_map = zeros(length(Yedges)-1,length(Xedges)-1);
        sptimes_unit_sec = spikeData.good_units.spikeTimes_usec{unit}/1e6; % spike times of the current unit in sec

        for ll = 1:n_land
            sptimes = sptimes_unit_sec( sptimes_unit_sec >= (t_land_sec(ll,1)) & sptimes_unit_sec <= (t_land_sec(ll,2)) ); % spike times 5 sec before and 5 sec after the current trial
            [~,idx] = arrayfun(@(sptime) min(abs(sptime-t)),sptimes); % Align spike times into the timestamps in behavioral data
            X = r(idx,1,bb)-r(idx,1,bbb); % x coordinates corresponding to each spike time
            Y = r(idx,2,bb)-r(idx,2,bbb); % y coordinates corresponding to each spike time
            [N,~,~] = histcounts2(Y,X,Yedges,Xedges); % row = Y, col = X
            spike_map(:,:) = spike_map(:,:) + N;
        end
        spike_map_flt = imgaussfilt(spike_map,flt_sigma);

        % Firing rate map
        rate_map(:,:,unit,bbb) = spike_map_flt ./ occupancy_map_flt(:,:,bbb);
    end
end

%% Save

save(fullfile(saveDir,'vector_coding.mat'),"Xedges",'Yedges','occupancy_map_flt','rate_map','bb')


