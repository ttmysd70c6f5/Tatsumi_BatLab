function Foraging_Task_v3(varargin)
%% Params and client init
global h_run h_stop h_bat h_feeders n_bats h_directory
global lightbar1 lightbar2 lightbar3 lightbar4
global feed1 feed2 feed3 feed4 ard trg h_d_rew
global h_feed1 h_feed2 h_feed3 h_feed4
global s_numbers;

ask_cdp = 1;
N = str2double(n_bats.String);

bhv_data.filename = h_directory.String;
bhv_data.fields = {'P1','P2','P3','P4','R1','R2','R3','R4','bat_ID,','activated_feeder','eventtime','outcome','rew_decr'};
bhv_data.trials = [];
save(bhv_data.filename,'bhv_data');
t_counter = 0;
start_time = tic;  %this is the time that each trial is referenced to

%Create port for getting CDP data from CDP_server (running on another
%instance of Matlab
cdp_client = tcpip('localhost', 30000, 'NetworkRole', 'client');
fopen(cdp_client);

%% Task body
while ~h_stop.Value
    rew_dur = cell2mat(h_feeders.Data(:,3));
    d_rew = str2double(h_d_rew.String);
    pf = cell2mat(h_feeders.Data(:,2));
    lb = ones(4,1);
    
    %Default state of the machine: check on beam breaks
    while all(lb) && ~h_stop.Value && ~h_feed1.Value && ~h_feed2.Value && ~h_feed3.Value && ~h_feed4.Value
        drawnow
        lb(1) = readDigitalPin(ard,lightbar1);
        lb(2) = readDigitalPin(ard,lightbar2);
        lb(3) = readDigitalPin(ard,lightbar3);
        lb(4) = readDigitalPin(ard,lightbar4);
        
        %check bat status from CDP server (feederID does not matter)
        for bat = 1:N
            fwrite(cdp_client,uint8([cell2mat(h_bat.Data(bat,2)), bat, 1, ask_cdp]),'uint8');
            y = fread(cdp_client,2,'uint8');    h_bat.Data{bat,2} = logical(y(1));
        end
        
    end
    
    %if one or more beam break occurs, check what's happening
    %disp('Event happened');
    if ~all(lb)
        disp(datetime('now'));
        interrupted_beams = find(~lb)                                                                  %check all the interrupted beams
        
        %check bat status from CDP server
        for bat = 1:N
            fwrite(cdp_client,uint8([cell2mat(h_bat.Data(bat,2)), bat, 1, ask_cdp]),'uint8');
            y = fread(cdp_client,2,'uint8');    h_bat.Data{bat,2} = logical(y(1));
        end
        
        %for all the interrupted beams
        for bb = 1:numel(interrupted_beams)
            
            active_feeder = interrupted_beams(bb);                                                      %putative active feeder
            fwrite(cdp_client,uint8([cell2mat(h_bat.Data(bat,2)), bat, active_feeder, ask_cdp]),'uint8');
            y = fread(cdp_client,2,'uint8');    who_triggered = y(2);                                   %who activated the feeder
            is_enabled = cell2mat(h_bat.Data(who_triggered,2));                                         %check if enabled to feed
            if is_enabled                                                                               %If enabled ->
                h_bat.Data{who_triggered,3} = cell2mat(h_bat.Data(who_triggered,3))+1;                  %Increase trial counter
                h_bat.Data{who_triggered,2} = false;                                                    %Disable bat that triggered feeder
                drawnow                                                                                 %Update gui
                %writeDigitalPin(ard,trg,1); pause(.05); writeDigitalPin(ard,trg,0);                     %Send 50 ms TTL from ARDUINO OUT
                if pf(active_feeder) >= (rand*100)                                                      %If probability allows
                    feeder_name = ['feed' num2str(active_feeder)];                                      %Get feeder name
                    start(eval(feeder_name));   pause(rew_dur(active_feeder));  stop(eval(feeder_name));%Feed
                    writeDigitalPin(ard,trg,1); pause(rew_dur(active_feeder));  writeDigitalPin(ard,trg,0); %Send TTL from ARDUINO OUT
                    h_bat.Data{who_triggered,4} = cell2mat(h_bat.Data(who_triggered,4))+1;              %Increase correct counter
                    h_feeders.Data{active_feeder,3} = max([cell2mat(h_feeders.Data(active_feeder,3))-d_rew,0.01]);  %Decrease reward (not more than 0.05)
                    drawnow                                                                             %update gui
                    
                    %Save behavioral file with 1 outcome
                    t_counter = t_counter + 1;  eventtime = toc(start_time);
                    bhv_data.trials(t_counter,:) = [pf', cell2mat(h_feeders.Data(:,3))', who_triggered, active_feeder, eventtime, 1, d_rew];
                    save(bhv_data.filename,'bhv_data');
                else
                    %Save behavioral file with 0 outcome
                    t_counter = t_counter + 1;  eventtime = toc(start_time);
                    bhv_data.trials(t_counter,:) = [pf', cell2mat(h_feeders.Data(:,3))', who_triggered, active_feeder, eventtime, 0, d_rew];
                    save(bhv_data.filename,'bhv_data');
                end
            end
        end
        
    end
    
    %If feed is requested
    drawnow
    if h_feed1.Value || h_feed2.Value || h_feed3.Value || h_feed4.Value
        feed;
    end
end

%% When stop is pressed
drawnow

%Kill CDP_server
ask_cdp = 0;
fwrite(cdp_client,uint8([1, 1, 1, ask_cdp]),'uint8');

%Terminate task
h_run.Value = 0;
h_stop.Value = 0;
disp('Stopped task');
end

