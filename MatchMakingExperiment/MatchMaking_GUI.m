function MatchMaking_GUI
%% Configure pins
ard=arduino('com5','Mega2560','BaudRate',115200);  %define arduino and port (yours might be com4 if I remember, it varies if you plug in to another port or change a computer

ttl_pin = 'A0'; %define the pin on arduino you want to use for ttl input
door_pin = 'D6';
OPEN_pin = 'D9';
CLOSE_pin = 'D12';

configurePin(ard,ttl_pin,'AnalogInput'); %configure the role of arduino pin
configurePin(ard,door_pin,'DigitalInput'); %configure the role of arduino pin
% configurePin(ard,OPEN_pin,'DigitalOutput'); %configure the role of arduino pin
configurePin(ard,CLOSE_pin,'DigitalOutput'); %configure the role of arduino pin

%% GUI
% h_fig = figure('name','4 Feeders task','ToolBar','none','DockControls','on','MenuBar','none');  
h_fig = uifigure('Name','4 Feeders task');
set(h_fig, 'units','normalized','Position',[0.2 0.1 0.5 0.6])

%hitting run starts the task
h_run = uicontrol(h_fig,'Style','togglebutton','String','Run','units','normalized',...
    'Position',[.05 .9 .4 .08],'Callback',@MatchMaking_Task,'fontsize',10,'fontweight','b','tag','run');

%stop task
h_stop = uicontrol(h_fig,'Style','togglebutton','String','stop','units','normalized',...
    'Position',[.55 .9 .4 .08],'fontsize',10,'fontweight','b','tag','stop');    

%next trial
h_next = uicontrol(h_fig,'Style','togglebutton','String','Next','units','normalized',...
    'Position',[.05 .8 .2 .08],'fontsize',10,'fontweight','b','tag','next');

%finish current trial
h_pause = uicontrol(h_fig,'Style','togglebutton','String','Pause','units','normalized',...
    'Position',[.30 .8 .2 .08],'fontsize',10,'fontweight','b','tag','next');

%Trial status
h_status = uicontrol(h_fig,'Style','radiobutton','String','Status','units','normalized',...
    'Position',[.55 .8 .1 .08],'fontsize',10,'fontweight','b','tag','next');

%default directory
h_dir = uicontrol(h_fig,'Style','text','String','Directory and File Name','units','normalized',...
   'Position',[.01 .74 .3 .05],'fontsize',9,'fontweight','b'); 
h_Dir = uicontrol(h_fig,'Style','edit','String',fullfile(cd,'Data'),'units','normalized',...
    'Position',[.01 .70 .65 .05],'fontsize',9,'fontweight','b');

h_name = uicontrol(h_fig,'Style','text','String','Rec Name','units','normalized',...
   'Position',[.67 .74 .3 .05],'fontsize',9,'fontweight','b'); 
h_Name = uicontrol(h_fig,'Style','edit','String',datestr(now,'mmddyy_HH_MM'),'units','normalized',...
    'Position',[.67 .70 .3 .05],'fontsize',9,'fontweight','b');

%Current session
h_trial = uicontrol(h_fig,'Style','text','String','Trial','units','normalized',...
   'Position',[.26 .64 .23 .05],'fontsize',9,'fontweight','b'); 
h_Trial = uicontrol(h_fig,'Style','edit','String','0','units','normalized',...
    'Position',[.26 .60 .23 .05],'fontsize',9,'fontweight','b');

% Current bat
h_ephys = uicontrol(h_fig,'Style','text','String','Ephys bat','units','normalized',...
   'Position',[.26 .54 .23 .05],'fontsize',9,'fontweight','b'); 
h_Ephys = uicontrol(h_fig,'Style','edit','String','','units','normalized',...
    'Position',[.26 .50 .23 .05],'fontsize',9,'fontweight','b');

% Matchmaking bat
h_bat = uicontrol(h_fig,'Style','text','String','Pair bat','units','normalized',...
   'Position',[.26 .44 .23 .05],'fontsize',9,'fontweight','b'); 
h_Bat = uicontrol(h_fig,'Style','edit','String','','units','normalized',...
    'Position',[.26 .40 .23 .05],'fontsize',9,'fontweight','b');

%% Task
function MatchMaking_Task(~,~)
    disp('Started task')

    % initialize a structure to save data
    data = [];
    trials = cell(0,5);
    t_counter = 0; % time counter
    trl_counter = 0; % trial counter
    savePath = fullfile(h_Dir.String,h_Name.String);

    % start counting timer
    start_time = tic;

    try
        while ~h_stop.Value
            ttl = readVoltage(ard,ttl_pin); %read the ttl voltage, should be 0 for baseline, and high voltage if there's apulse
            
            % Get the current door status
            door = readDigitalPin(ard,door_pin); % Normally open (OPEN = 1)
    
            % Update LED
            if door
                % writeDigitalPin(ard,OPEN_pin,1)
                writeDigitalPin(ard,CLOSE_pin,0)
            else
                % writeDigitalPin(ard,OPEN_pin,0)
                writeDigitalPin(ard,CLOSE_pin,1)
            end
    
            % save ttl voltage and door status
            save_signal()

            if h_next.Value == 1
                pre_trial = str2double(h_Trial.String);
                h_Trial.String = num2str(pre_trial+1); % update trial
                save_trial('start')
                h_next.Value = 0;
                h_status.Value = 1;
            end

            if h_pause.Value == 1
                save_trial('stop')
                h_pause.Value = 0;
                h_status.Value = 0;
            end
        end

        % Kill task
        stopTask()
    catch
        % Kill task
        stopTask()
    end

    function save_signal
        t_counter = t_counter + 1;
        eventtime = toc(start_time);
        data(t_counter,:) = [eventtime,door,ttl];
    end

    function save_trial(action)
        trl_counter = trl_counter + 1;

        eventtime = toc(start_time);
        trial = str2double(h_Trial.String);
        ephys_bat = h_Ephys.String;
        pair_bat = h_Bat.String;
        trials = [trials; {eventtime,trial,action,ephys_bat,pair_bat}];
    end

    function stopTask
        % stop
        drawnow
        
        %Terminate task
        h_run.Value = 0;
        h_next.Value = 0;
    
        writeDigitalPin(ard,OPEN_pin,0)
        writeDigitalPin(ard,CLOSE_pin,0)
    
        save(savePath,'data','trials')
    
        disp('Stopped task');

        clear all
    end
end

end