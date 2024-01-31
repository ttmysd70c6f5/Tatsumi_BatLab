function record_audio_AF(varargin)

batName = [];
dateSesh = datestr(now,'yymmdd');

% User inputs overrides
nparams=length(varargin);
for i=1:2:nparams
    switch lower(varargin{i})
        case 'batname'
            batName=varargin{i+1};
        case 'datesesh'
            dateSesh=varargin{i+1};
    end
end
%records audio and sync pulses
global h_directory_aud h_run_aud h_stop_aud
global IDsound fs rec_dur input_channels inter_file_time
global homeDir batName dateSesh

%create gui
h_fig = figure('name','Continuous Audio Recording','Position',[65 143 300 500],'ToolBar','none','DockControls','on','MenuBar','none');  

%default directory
homeDir = ['D:\aforli\Audio_Data' batName '\' dateSesh];
h_dir = uicontrol(h_fig,'Style','text','String','Directory','units','normalized',...
   'Position',[.01 .84 .3 .05],'fontsize',9,'fontweight','b'); 
h_directory_aud = uicontrol(h_fig,'Style','edit','String',homeDir,'units','normalized',...
    'Position',[.01 .80 .96 .05],'fontsize',9,'fontweight','b');

%hitting run starts the task
h_run_aud = uicontrol(h_fig,'Style','togglebutton','String','Run','units','normalized',...
    'Position',[.01 .9 .3 .1],'Callback',@run_rec_audio_AF,'fontsize',10,'fontweight','b','tag','run');

%stop task
h_stop_aud = uicontrol(h_fig,'Style','togglebutton','String','Stop','units','normalized',...
    'Position',[.32 .9 .3 .1],'fontsize',10,'fontweight','b','tag','run');


%audio
IDsound = 'ASIO HDSPe FX';
fs = 192000;
rec_dur = 12; %set to 30 seconds (trials will be stored every 30)
input_channels = 0:7; %0-5 are the microphones and 6-7 TTL inputs (only collect mics on feeders 3 and 4 + side mics)

inter_file_time = 10; %this is the time between saving each file
%%above is all that is necessary for recording sound
