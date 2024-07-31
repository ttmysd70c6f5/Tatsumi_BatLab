function Foraging_Task(varargin)
%% Params and client init
global h_directory h_run h_pause h_stop n_bats h_bat h_feeders h_session h_d_rew h_Name h_Exp h_Alt h_Alt_1 h_Alt_2 h_Alt_now
global lightbar1 lightbar2 lightbar3 lightbar4
global feed1 feed2 feed3 feed4 ard trg
global h_feed1 h_feed2 h_feed3 h_feed4 
global h_Rew1 h_Rew2 h_Rew3 h_Rew4
global s_numbers bat_ids

ask_cdp = 1;
N = str2double(n_bats.String); % # of bats
n_fd = 4; % # of feeders

bhv_data.filename = h_directory.String; % directory to save data
bhv_data.fields = {'P1','P2','P3','P4','R1','R2','R3','R4','S1','S2','S3','S4','bat_ID','activated_feeder','eventtime','outcome','rew_decr','session','sw_th'};
bhv_data.trials = [];
% bhv_data.bats = bat_ids;
bhv_data.bats = h_bat.Data(:,2);
bhv_data.s_numbers = s_numbers;
bhv_data.tags = h_bat.Data(:,1);
bhv_data.recname = h_Name.String; % recording name
bhv_data.recdate = datestr(now,'mmddyy_HHMM'); % rec date
bhv_data.reward = {h_Rew1.String; h_Rew2.String; h_Rew3.String; h_Rew4.String};
bhv_data.ephys = bhv_data.bats(cell2mat(h_bat.Data(:,6))); % ephys bats
bhv_data.region = h_bat.Data(cell2mat(h_bat.Data(:,6)),7); % ephys region
bhv_data.exptype = h_Exp.String(h_Exp.Value); % experiment type
bhv_data.switch = true;
sw_1 = str2double(h_Alt_1.String);
sw_2 = str2double(h_Alt_2.String);
bhv_data.switch_fd = [sw_1;sw_2];
save(bhv_data.filename,'bhv_data');
t_counter = 0;
start_time = tic;  %this is the time that each trial is referenced to
pause_state = false; % true if a task is pausing
h_feeders.Data(:,3) = h_feeders.Data(:,2); % Initialize feeder states

%Create port for getting CDP data from CDP_server (running on another
%instance of Matlab
cdp_client = tcpip('localhost', 30000, 'NetworkRole', 'client');
fopen(cdp_client);

%% Task body
while ~h_stop.Value
    lb = ones(4,1); % state of lightbar (0=beam breaked) [4,1]
    fd_enabled = cell2mat(h_feeders.Data(:,2)); % Only enabled feeder is considered        

    %Default state of the machine: check on beam breaks
    while all(lb(fd_enabled)) && ~h_pause.Value && ~h_stop.Value && ~h_feed1.Value && ~h_feed2.Value && ~h_feed3.Value && ~h_feed4.Value
        % Update settings
        rew_dur = cell2mat(h_feeders.Data(:,5)); % Reward duration [4,1]
        d_rew = str2double(h_d_rew.String); % Reward decrement [1]
        pf = cell2mat(h_feeders.Data(:,4)); % Reward probability [4,1]
        fd_enabled = cell2mat(h_feeders.Data(:,2)); % Only enabled feeder is considered        
        fd_state = cell2mat(h_feeders.Data(:,3)); % 1 if the feeder is on
        session = str2num(h_session.String); % current session
        
        drawnow % update GUI

        if fd_enabled(1);   lb(1) = readDigitalPin(ard,lightbar1);  end
        if fd_enabled(2);   lb(2) = readDigitalPin(ard,lightbar2);  end
        if fd_enabled(3);   lb(3) = readDigitalPin(ard,lightbar3);  end
        if fd_enabled(4);   lb(4) = readDigitalPin(ard,lightbar4);  end

        %check bat status from CDP server (feederID does not matter)
        for bat = 1:N
            fwrite(cdp_client,uint8([1, cell2mat(h_bat.Data(bat,3)), bat, 1, ask_cdp]),'uint8'); % type, state, bat, feeder, ask_cdp
            y = fread(cdp_client,1,'uint8');    h_bat.Data{bat,3} = logical(y);
        end

        % Check feeder status from CDP server
        for fd = 1:n_fd
            if fd_enabled(fd) && ~fd_state(fd)
                fwrite(cdp_client,uint8([3, fd_state(fd), 1, fd, ask_cdp]),'uint8'); % type, state, bat, feeder, ask_cdp
                y = fread(cdp_client,1,'uint8');    h_feeders.Data{fd,3} = logical(y);
%                 disp(sprintf('%d: Enabled feeder%d',y,fd)) 
            end
        end
    end


    %if one or more beam break occurs, check what's happening
    %disp('Event happened');
    if ~all(lb(fd_enabled)) && ~h_pause.Value && ~h_stop.Value
        interrupted_beams = find(~lb);                                                                  %check all the interrupted beams
        fprintf('%s: feeder %d is blocked\n',datetime('now'),interrupted_beams);

        %check bat status from CDP server
        for bat = 1:N
            fwrite(cdp_client,uint8([1, cell2mat(h_bat.Data(bat,3)), bat, 1, ask_cdp]),'uint8');
            y = fread(cdp_client,1,'uint8');    h_bat.Data{bat,3} = logical(y); % if enabled to feed
        end
        
        %for all the interrupted beams
        for bb = 1:numel(interrupted_beams)
            active_feeder = interrupted_beams(bb);                                                          %putative active feeder
          
            % Check who activated the feeder
            fwrite(cdp_client,uint8([2, true, 1, active_feeder, ask_cdp]),'uint8'); % Only feeder ID does matter
            y = fread(cdp_client,1,'uint8');    who_triggered = y;                                   %who activated the feeder

            is_enabled = cell2mat(h_bat.Data(who_triggered,3));                                         %check if enabled to feed
            if is_enabled                                                                               %If enabled ->
                h_bat.Data{who_triggered,4} = cell2mat(h_bat.Data(who_triggered,4))+1;                  %Increase trial counter
                h_bat.Data{who_triggered,3} = false;                                                    %Disable bat that triggered feeder
                drawnow                                                                                 %Update gui
                %writeDigitalPin(ard,trg,1); pause(.05); writeDigitalPin(ard,trg,0);                     %Send 50 ms TTL from ARDUINO OUT
                if fd_state(active_feeder) && fd_enabled(active_feeder)                                       %if active feeder is on
                    h_feeders.Data{active_feeder,3} = false;                                                    %Turn off the active feeder
                    
                    if pf(active_feeder) >= (rand*100)                                                      %If probability allows
                        h_Alt_now.String = num2str(str2double(h_Alt_now.String)+1);                         % Count feeder activation
                        
                        feeder_name = ['feed' num2str(active_feeder)];                                      %Get feeder name
                        start(eval(feeder_name));   pause(rew_dur(active_feeder));  stop(eval(feeder_name));%Feed
                        writeDigitalPin(ard,trg,1); pause(rew_dur(active_feeder));  writeDigitalPin(ard,trg,0); %Send TTL from ARDUINO OUT
                        h_bat.Data{who_triggered,5} = cell2mat(h_bat.Data(who_triggered,5))+1;              %Increase correct counter
                        h_feeders.Data{active_feeder,5} = max([cell2mat(h_feeders.Data(active_feeder,5))-d_rew,0.01]);  %Decrease reward (not less than 0.01)
                        drawnow                                                                             %update gui
                        
                        %Save behavioral file with 1 outcome
                        t_counter = t_counter + 1;  eventtime = toc(start_time);
                        bhv_data.trials(t_counter,:) = [pf', cell2mat(h_feeders.Data(:,5))',cell2mat(h_feeders.Data(:,6))', who_triggered, active_feeder, eventtime, 1, d_rew, session, str2double(h_Alt.String)];
                        save(bhv_data.filename,'bhv_data');
                    else
                        %Save behavioral file with 0 outcome
                        t_counter = t_counter + 1;  eventtime = toc(start_time);
                        bhv_data.trials(t_counter,:) = [pf', cell2mat(h_feeders.Data(:,5))', cell2mat(h_feeders.Data(:,6))', who_triggered, active_feeder, eventtime, 0, d_rew, session, str2double(h_Alt.String)];
                        save(bhv_data.filename,'bhv_data');
                    end
                end
            end
        end
        
        if str2double(h_Alt.String) ~= 0 && str2double(h_Alt_now.String) >= str2double(h_Alt.String)                        % If reached threshold
            pr_1_now = h_feeders.Data(sw_1,4);
            pr_2_now = h_feeders.Data(sw_2,4);
            
            h_feeders.Data(sw_1,4) = pr_2_now; % switch probabilities between feeders
            h_feeders.Data(sw_2,4) = pr_1_now;
            
            h_Alt_now.String = '0'; % reset counter
        end
    end

    %If feed is requested
    drawnow
    if h_feed1.Value || h_feed2.Value || h_feed3.Value || h_feed4.Value
        % feed;
        feed_button();
    end


    % Pause task: Reset bat and lightbar status and record the time
    if h_pause.Value && ~h_stop.Value
        for bat = 1:N
            h_bat.Data{bat,3} = true; % reset bat status
        end
        lb = ones(4,1); % reset lightbar status

        %Save behavioral file with 0 for bat and feeder
        if ~pause_state
            t_counter = t_counter + 1;  eventtime = toc(start_time);
            bhv_data.trials(t_counter,:) = [pf', cell2mat(h_feeders.Data(:,5))',cell2mat(h_feeders.Data(:,6))', -1, 0, eventtime, 0, d_rew, session, str2double(h_Alt.String)]; % -1 in the 9th column when task paused
            save(bhv_data.filename,'bhv_data');
            pause_state = true; % Record pause only once
        end
    end


    % Resume task: record the time
    if ~h_pause.Value && ~h_stop.Value && pause_state
        %Save behavioral file with 0 for bat and feeder
        t_counter = t_counter + 1;  eventtime = toc(start_time);
        session = str2num(h_session.String); % current session
        bhv_data.trials(t_counter,:) = [pf', cell2mat(h_feeders.Data(:,5))', cell2mat(h_feeders.Data(:,6))', -2, 0, eventtime, 0, d_rew, session, str2double(h_Alt.String)]; % -2 in the 9th column when task resumed
        save(bhv_data.filename,'bhv_data');
        pause_state = false; % Enable recording pause for next time
    end
end

%% When stop is pressed
drawnow

%Kill CDP_server
ask_cdp = 0;
fwrite(cdp_client,uint8([1, 1, 1, 1, ask_cdp]),'uint8');

%Terminate task
h_run.Value = 0;
h_pause.Value = 0;
h_stop.Value = 0;
disp('Stopped task');
end

