function exportCEBRA(BhvData, spikeData, mtData, bname_ephys, saveFile)
% Export neural and behavioral data to python for running CEBRA

%%% Initialize the struct to output
out = struct([]);

%%% Extract data only for the bat that the neural activities were recorded
bb = find(ismember(mtData.bname,bname_ephys));

%%% Load behavioral data
bflying = BhvData.bflying; % flying state
t = BhvData.t; % time in CDP
n_bats = length(mtData.bname);
R = zeros(length(t),n_bats);
for bbb = 1:length(mtData.bname)
    [~, R(:,bbb), ~] = curvature(BhvData.r(:,:,bbb));
end
BhvData.R = R;

%%% Flight detection
flight_start_idx    =   find(diff([bflying(1,bb);bflying(:,bb)])==1); % The index of the variable t for the start of flights
flight_end_idx      =   find(diff([bflying(:,bb);bflying(end,bb)])==-1); % The index of the variable t for the end of flights
% Omit the initial/last flight if the recording was started/stopped during a flight
if bflying(1,bb)~=1 % if the bat was not flying when the recording started
    if bflying(end,bb)~=1 % if the bat was not flying when the recording finished
        % no correction
    else % if the bat was landing when the recording finished
        flight_start_idx(end) = []; % omit the final takeoff
    end
elseif bflying(1,bb)==1 % if the bat was flying when the recording started
    if bflying(end,bb)~=1 % if the bat was not flying when the recording finished
        flight_end_idx(1) = []; % omit the initial landing
    else % if the bat was flying when the recording finished
        flight_start_idx(end) = []; % omit the final takeoff
        flight_end_idx(1) = []; % omit the initial landing
    end
end

flight_start_sec = t(flight_start_idx); % timestamps of the starts of flights in sec
flight_end_sec = t(flight_end_idx); % timestamps of the ends of flights in sec

n_flight = length(flight_start_idx); % number of flight
n_unit = spikeData.numGoodUnits; % number of good units

%%% FLY: Extract spike times during the current flight
sptimes_sec = [];
sp_unit = [];
bhv_times_sec = [];
position = [];
velocity = [];
accel = [];
curv = [];
for flight = 2:n_flight-1  
    for unit = 1:n_unit % the current unit
        % spike times of the current unit during the current flight (sec)
        sptimes_unit_sec = spikeData.good_units.spikeTimes_usec{unit}/1e6;
        sptimes_unit_idx = sptimes_unit_sec >= (flight_start_sec(flight)) & sptimes_unit_sec <= (flight_end_sec(flight));
        sptimes_unit_sec = sptimes_unit_sec(sptimes_unit_idx);
        sptimes_sec = [sptimes_sec; sptimes_unit_sec];
        sp_unit = [sp_unit; repmat(unit,length(sptimes_unit_sec),1)];
    end
    bhv_idx = flight_start_idx(flight):flight_end_idx(flight);
    bhv_times_sec = [bhv_times_sec; t(bhv_idx)];
    position = [position; BhvData.r(bhv_idx,:,bb)];
    velocity = [velocity; BhvData.v_abs(bhv_idx,bb)];
    accel = [accel; BhvData.ac_abs(bhv_idx,bb)];
    curv = [curv; BhvData.R(bhv_idx,bb)]; % curvature
    % curv = R;
end
% save
out(1).Neural(1).Fly(1).SpikeTimes_sec = sptimes_sec;
out(1).Neural(1).Fly(1).UnitIDs = sp_unit;
out(1).Behavior(1).Fly(1).BhvTimes_sec = bhv_times_sec;
out(1).Behavior(1).Fly(1).position = position;
out(1).Behavior(1).Fly(1).velocity = velocity;
out(1).Behavior(1).Fly(1).accel = accel;
out(1).Behavior(1).Fly(1).curvature = curv;
    

%%% PRE: Extract spike times at the resting position before the current flight
sptimes_sec = [];
sp_unit = [];
bhv_times_sec = [];
position = [];
velocity = [];
accel = [];
curv = [];
for flight = 2:n_flight-1  
    for unit = 1:n_unit % the current unit
        % spike times of the current unit during the current flight (sec)
        sptimes_unit_sec = spikeData.good_units.spikeTimes_usec{unit}/1e6;
        if (flight_start_sec(flight)) - (flight_end_sec(flight-1)) <= 30 % if the previous landing is less than 30 sec
            sptimes_unit_idx = sptimes_unit_sec >= (flight_end_sec(flight-1)) & sptimes_unit_sec <= (flight_start_sec(flight));
            sptimes_unit_sec = sptimes_unit_sec(sptimes_unit_idx);
        else
            sptimes_unit_idx = sptimes_unit_sec >= (flight_start_sec(flight)-30) & sptimes_unit_sec <= (flight_start_sec(flight));
            sptimes_unit_sec = sptimes_unit_sec(sptimes_unit_idx);
        end
        sptimes_sec = [sptimes_sec; sptimes_unit_sec];
        sp_unit = [sp_unit; repmat(unit,length(sptimes_unit_sec),1)];
    end
    if (flight_start_sec(flight)) - (flight_end_sec(flight-1)) <= 30 % if the previous landing is less than 30 sec
        bhv_idx = flight_end_idx(flight-1):flight_start_idx(flight);
    else
        bhv_idx = flight_start_idx(flight)-30:flight_start_idx(flight);
    end
    bhv_times_sec = [bhv_times_sec; t(bhv_idx)];
    position = [position; BhvData.r(bhv_idx,:,bb)];
    velocity = [velocity; BhvData.v_abs(bhv_idx,bb)];
    accel = [accel; BhvData.ac_abs(bhv_idx,bb)];
    curv = [curv; BhvData.R(bhv_idx,bb)]; % curvature
end
% save
out(1).Neural(1).Pre(1).SpikeTimes_sec = sptimes_sec;
out(1).Neural(1).Pre(1).UnitIDs = sp_unit;
out(1).Behavior(1).Pre(1).BhvTimes_sec = bhv_times_sec;
out(1).Behavior(1).Pre(1).position = position;
out(1).Behavior(1).Pre(1).velocity = velocity;
out(1).Behavior(1).Pre(1).accel = accel;
out(1).Behavior(1).Pre(1).curvature = curv;


%%% POST: Extract spike times at the resting position after the current flight
sptimes_sec = [];
sp_unit = [];
bhv_times_sec = [];
position = [];
velocity = [];
accel = [];
curv = [];
for flight = 2:n_flight-1  
    for unit = 1:n_unit % the current unit
        % spike times of the current unit during the current flight (sec)
        sptimes_unit_sec = spikeData.good_units.spikeTimes_usec{unit}/1e6;
        if (flight_start_sec(flight+1)) - (flight_end_sec(flight)) <= 30 % if the previous landing is less than 30 sec
            sptimes_unit_idx = sptimes_unit_sec >= (flight_end_sec(flight)) & sptimes_unit_sec <= (flight_start_sec(flight+1));
            sptimes_unit_sec = sptimes_unit_sec(sptimes_unit_idx);
        else
            sptimes_unit_idx = sptimes_unit_sec >= (flight_end_sec(flight)) & sptimes_unit_sec <= (flight_end_sec(flight)+30);
            sptimes_unit_sec = sptimes_unit_sec(sptimes_unit_idx);
        end
        sptimes_sec = [sptimes_sec; sptimes_unit_sec];
        sp_unit = [sp_unit; repmat(unit,length(sptimes_unit_sec),1)];
    end
    if (flight_start_sec(flight+1)) - (flight_end_sec(flight)) <= 30 % if the previous landing is less than 30 sec
        bhv_idx = flight_end_idx(flight):flight_start_idx(flight+1);
    else
        bhv_idx = flight_end_idx(flight):flight_end_idx(flight)+30;
    end
    bhv_times_sec = [bhv_times_sec; t(bhv_idx)];
    position = [position; BhvData.r(bhv_idx,:,bb)];
    velocity = [velocity; BhvData.v_abs(bhv_idx,bb)];
    accel = [accel; BhvData.ac_abs(bhv_idx,bb)];
    curv = [curv; BhvData.R(bhv_idx,bb)]; % curvature
end
% save
out(1).Neural(1).Post(1).SpikeTimes_sec = sptimes_sec;
out(1).Neural(1).Post(1).UnitIDs = sp_unit;
out(1).Behavior(1).Post(1).BhvTimes_sec = bhv_times_sec;
out(1).Behavior(1).Post(1).position = position;
out(1).Behavior(1).Post(1).velocity = velocity;
out(1).Behavior(1).Post(1).accel = accel;
out(1).Behavior(1).Post(1).curvature = curv;


% saveFile = fullfile(rootdir,'analysis',sprintf('cebra_%s_%s_%s.mat',bname_ephys,trjName,mtData.recName));
save(saveFile,'out')