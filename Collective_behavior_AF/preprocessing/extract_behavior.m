function extract_behavior(prm)
sessionFile = dir(fullfile(pwd,'*/*cdp*'));

%% Parameters
fl_offset = prm.fl_offset;                                           % Consider flights after this value (secs)


%% Process for all sessions
fl_cnt = 1; % Flight count for assigining shared ID
bout_cnt = 1; % Bout count for assigining shared ID

for ss = 1:length(sessionFile)
    sessionDir = sessionFile(ss).folder;

    %% Load data
    %=== Behavioral data and related
    BHV1_file = dir(fullfile(sessionDir, 'Ext_Beh*','Extracted_Behavior*'));                              % Preprocessed Behavioral Data
    load(fullfile(BHV1_file.folder,BHV1_file.name),'a','a_abs','angle','bat_clr','bat_nms','bat_pair_nms','bat_pairs',...
        'batdate','bflying','Fs','f_num','Group_name','n_tags','r','r_lim','t','T','v','v_abs','v_th','wBeats');        % Processed Behavioral Data
    
    BHV2_file = dir(fullfile(sessionDir, 'Ext_Beh*','Analysis*/Analyzed_Behavior*'));
    load(fullfile(BHV2_file.folder,BHV2_file.name));        % Analyzed behavior
    
    bout_file = dir(fullfile(sessionDir, 'Ext_Beh*','bouts*'));       load(fullfile(bout_file.folder,bout_file.name));        % flight bout
    
    bpos_file = dir(fullfile(sessionDir, 'Ext_Beh*','bpos*'));
    load(fullfile(bpos_file.folder,bpos_file.name));
    
    % Extract session number and date
    dirStr = strsplit(sessionDir,'\');
    DatasetNum = dirStr{end-1};        RecDate = dirStr{end};
    
    %% Summary of behaviors
    % Extract_behavior extracts flights and bouts from the session specified
    % by sessionDir.
    
    % Extract flights
    bat_id          =       (1:n_tags)';
    dataset_num     =       repmat(DatasetNum,[n_tags,1]);
    rec_date        =       repmat(RecDate,[n_tags,1]);           
    fl_num          =       zeros(n_tags,1);
    
    for nb = 1:n_tags
        [takeoff, landing] = findchanges(bflying(:,nb),0,1);                    % Extract takeoff/landing timings
        t_toff = t(takeoff);                                                 % Time of takeoff
        t_land = t(landing);                                                 % Time of landing   
        if length(t_toff) > length(t_land)
            t_toff = t_toff(1:end-1);                                     % Remove the final takeoff in the case the session interupts the final flight
        elseif length(t_toff) < length(t_land)
            error('Mismatch of # of takeoffs and landings')
        end
        t_toff = t_toff(t_toff > fl_offset);                           % Exclude the initial flights
    
        fl_num(nb) = length(t_toff);                                         % Number of flights
    end
    
    fl_rate         =       fl_num / T * Fs;                                        % Flight rate in session
    
    % Create table and save
    behaviors = table(bat_id, dataset_num, rec_date, fl_num, fl_rate);

    %% Extract flights
    t_toff      =       [];             % Time of takeoff
    t_land      =       [];             % Time of landing
    fl_len       =       [];             % Length of flight
    ts_toff     =       [];             % Timestamp of takeoff
    ts_land     =       [];             % Timestamp of landing
    bat_id      =       [];             % Bat number
    r_toff      =       [];             % location when takeoff
    r_land      =       [];             % location when landing
    pos_toff    =       [];             % position cluster when takeoff
    pos_land    =       [];         % position cluster when landing
    
    for nb = 1:n_tags
        [takeoff, landing] = findchanges(bflying(:,nb),0,1);                    % Extract takeoff/landing timings
        t_toff_now       =       t(takeoff);                       
        t_land_now       =       t(landing);                        
        ts_toff_now      =       find(takeoff);                    
        ts_land_now      =       find(landing);                    
        if length(t_toff_now) > length(t_land_now)
            t_toff_now      =       t_toff_now(1:end-1);                        % Remove the final takeoff in the case the session interupts the final flight
            ts_toff_now     =       ts_toff_now(1:end-1);
        elseif length(t_toff) < length(t_land)
            error('Mismatch of # of takeoffs and landings')
        end
    
        t_toff          =       [t_toff; t_toff_now];
        t_land          =       [t_land; t_land_now];
        fl_len           =       [fl_len; t_land_now - t_toff_now];
        ts_toff         =       [ts_toff; ts_toff_now];
        ts_land         =       [ts_land; ts_land_now];
        r_toff          =       [r_toff; r(ts_toff_now,:,nb)];
        r_land          =       [r_land; r(ts_land_now,:,nb)];
        pos_toff        =       [pos_toff; bpos(ts_toff_now-1,nb)];
        pos_land        =       [pos_land; bpos(ts_land_now+1,nb)];
        bat_id          =       [bat_id; repmat(nb, [sum(takeoff),1])];
    
    end
    
    dataset_num     =       repmat(DatasetNum,[size(t_toff,1),1]);              % Dataset number
    rec_date        =       repmat(RecDate,[size(t_toff,1),1]);                 % Session date

    % Create table
    flights       =       table(dataset_num, rec_date, bat_id, t_toff, t_land, ts_toff, ts_land, fl_len,...
        r_toff, r_land, pos_toff, pos_land);
    flights       =       flights(flights.t_toff > fl_offset, :);         % Exclude initial takeoffs

    % Assign global id
    flid            =       [fl_cnt:fl_cnt-1+size(flights,1)]';                                % Shared flight ID
    flights.flid    =       flid;
    flights = movevars(flights, 'flid', 'Before', 'dataset_num');
    fl_cnt          =       fl_cnt+size(flights,1);                              % Start ID for next session

    %% Extract all bouts
    t_toff          =       [];             % Time of takeoff
    t_land          =       [];             % Time of landing
    bout_len          =       [];             % Length of flight
    ts_toff         =       [];             % Timestamp of takeoff
    ts_land         =       [];             % Timestamp of landing
    bat_id          =       [];             % Bat number
    r_toff          =       [];             % location when takeoff
    r_land          =       [];             % location when landing
    pos_toff        =       [];             % position cluster when takeoff
    pos_land        =       [];             % position cluster when landing
    num_fl_bout     =       [];             % Number of flights during bouts
    rate_fl_bout    =       [];             % Flight rate during bouts
    
    for nb = 1:n_tags
        [takeoff, landing] = findchanges(bout(:,nb),0,1);                    % Extract the first takeoff and the last landing in bouts
        t_toff_now       =       t(takeoff);                       
        t_land_now       =       t(landing);                        
        ts_toff_now      =       find(takeoff);                    
        ts_land_now      =       find(landing);                    
        if length(t_toff_now) > length(t_land_now)
            t_toff_now      =       t_toff_now(1:end-1);                        % Remove the final takeoff in the case the session interupts the final flight
            ts_toff_now     =       ts_toff_now(1:end-1);
        elseif length(t_toff) < length(t_land)
            error('Mismatch of # of takeoffs and landings')
        end
    
        [takeoff_fl, ~] = findchanges(bflying(:,nb),0,1);                   % takeoff timings of all flights
        for nf = 1:length(t_toff_now)
            cnt_fl       = sum(takeoff_fl(ts_toff_now(nf):ts_land_now(nf)));    % Count the number of flights during a bout
            num_fl_bout  = [num_fl_bout; cnt_fl];
            rate_fl_bout = [rate_fl_bout; cnt_fl ./ (t_land_now(nf)-t_toff_now(nf))];     % Compute the rate of flights during a bout
        end
    
        t_toff          =       [t_toff; t_toff_now];
        t_land          =       [t_land; t_land_now];
        ts_toff         =       [ts_toff; ts_toff_now];
        ts_land         =       [ts_land; ts_land_now];
        r_toff          =       [r_toff; r(ts_toff_now,:,nb)];
        r_land          =       [r_land; r(ts_land_now,:,nb)];
        pos_toff        =       [pos_toff; bpos(ts_toff_now-1,nb)];
        pos_land        =       [pos_land; bpos(ts_land_now+1,nb)];
        bat_id          =       [bat_id; repmat(nb, [sum(takeoff),1])];
    
    end
    
    bout_len           =       t_land - t_toff;
    dataset_num     =       repmat(DatasetNum,[size(t_toff,1),1]);              % Dataset number
    rec_date        =       repmat(RecDate,[size(t_toff,1),1]);                 % Session date

    % Create table and save
    bouts           =       table(dataset_num, rec_date, bat_id, t_toff, t_land, ts_toff, ts_land, bout_len,...
        num_fl_bout, rate_fl_bout, r_toff, r_land, pos_toff, pos_land);
    bouts           =       bouts(bouts.t_toff > fl_offset, :);         % Exclude initial takeoffs

    % Assign shared id
    boutid          =       [bout_cnt:bout_cnt-1+size(bouts,1)]';                                % Shared bout ID
    bouts.boutid    =       boutid;
    bouts = movevars(bouts, 'boutid', 'Before', 'dataset_num');
    bout_cnt        =       bout_cnt+size(bouts,1);                        % start Id for next session

    %% Save
    mkdir(fullfile(sessionDir,'Analysis_TY'))
    save(fullfile(sessionDir,'Analysis_TY','Behaviors_Extracted.mat'),'behaviors','flights','bouts')
end