function Foraging_Task_GUI
clc;
clear all;
clear global;

%Load configuration file
load('FR_config.mat');

%Declare global variables
global h_directory h_run h_pause h_stop n_bats h_bat h_feeders h_session h_Name h_Exp
global lightbar1 lightbar2 lightbar3 lightbar4
global feed1 feed2 feed3 feed4 ard trg
global rewardspeed h_d_rew
global h_feed1 h_feed2 h_feed3 h_feed4 h_reload
global s_numbers bat_ids;
global h_Rew1 h_Rew2 h_Rew3 h_Rew4
% global fed_override;

%Insert here the serial numbers of the tags (ORDERED)
% s_numbers = load('Serial_numbers_in_use.txt');
s_table = readtable('serial_numbers.csv','Format','%d %s %s'); s_table = s_table(~strcmp(s_table.Bat_id,''),:);
bat_ids = s_table.Bat_id;
s_numbers = double(s_table.Serial_number);
tag_names = s_table.Name;
% bat_names = ['Bt00'; 'Bt01'; 'Bt02'; 'Bt03'; 'Bt04'; 'Bt05'; 'Bt06'; 'Bt07'; 'Bt08'; 'Bt09'; 'Bt10'; 'Bt11'; 'Bt12'; 'Bt13'; 'Bt14'; 'Bt15'];
% bat_names = bat_names(1:length(s_numbers),:);
rewardspeed = 1; 
session_init = 1; % initialize session

%Initialize and configure arduino
ard=arduino('com3','uno','Libraries','Adafruit\MotorShieldV2'); % create an arduino object
shield=addon(ard,'Adafruit\MotorShieldV2'); % creates an add-on connection to Adafruit Motor Shield V2 connected to the Arduino hardware
trg = 'D7';         configurePin(ard,trg,'pullup');     configurePin(ard,trg,'DigitalOutput');      writeDigitalPin(ard,trg,0); % the pin D7 -> 'pullup' and 'DigitalOutput' mode, value 0 
lightbar1 = 'D3';   configurePin(ard,lightbar1,'pullup'); % Set D3 to 'pullup' mode
lightbar2 = 'D2';   configurePin(ard,lightbar2,'pullup'); % set D2 to 'pullup' mode
lightbar3 = 'D4';   configurePin(ard,lightbar3,'pullup'); % set D4 to 'pullup' mode
lightbar4 = 'D5';   configurePin(ard,lightbar4,'pullup'); % set D5 to 'pullup' mode
MOT1 = 2;           feed1=dcmotor(shield,MOT1);     feed1.Speed = rewardspeed; % configure the feeder 2
MOT2 = 1;           feed2=dcmotor(shield,MOT2);     feed2.Speed = rewardspeed; % configure the feeder 1
MOT3 = 3;           feed3=dcmotor(shield,MOT3);     feed3.Speed = rewardspeed; % configure the feeder 3
MOT4 = 4;           feed4=dcmotor(shield,MOT4);     feed4.Speed = rewardspeed; % configure the feeder 4
disp('--ARDUINO correctly initialized--');

%Set seed for random number generation
rand('state',sum(100*clock));

%% Create GUI
h_fig = figure('name','4 Feeders task','ToolBar','none','DockControls','on','MenuBar','none');  
% h_fig = uifigure('Name','4 Feeders task');
set(h_fig, 'units','normalized','Position',[0.05 0.3 0.5 0.5])
% [65 143 700 500]

%hitting run starts the task
h_run = uicontrol(h_fig,'Style','togglebutton','String','Run','units','normalized',...
    'Position',[.01 .9 .3 .1],'Callback',@Foraging_Task,'fontsize',10,'fontweight','b','tag','run');

%pause task
h_pause = uicontrol(h_fig,'Style','togglebutton','String','Pause','units','normalized',...
    'Position',[.32 .9 .3 .1],'fontsize',10,'fontweight','b','tag','run');

%stop task
h_stop = uicontrol(h_fig,'Style','togglebutton','String','Stop','units','normalized',...
    'Position',[.63 .9 .3 .1],'fontsize',10,'fontweight','b','tag','run');


%default directory
h_dir = uicontrol(h_fig,'Style','text','String','Directory and File Name','units','normalized',...
   'Position',[.01 .84 .3 .05],'fontsize',9,'fontweight','b'); 
h_directory = uicontrol(h_fig,'Style','edit','String',fullfile(cd,'Data', datestr(now,'mmddyy_HH_MM')),'units','normalized',...
    'Position',[.01 .80 .65 .05],'fontsize',9,'fontweight','b');

h_name = uicontrol(h_fig,'Style','text','String','Rec Name','units','normalized',...
   'Position',[.67 .84 .3 .05],'fontsize',9,'fontweight','b'); 
h_Name = uicontrol(h_fig,'Style','edit','String','','units','normalized',...
    'Position',[.67 .80 .3 .05],'fontsize',9,'fontweight','b');

%Number of bats
h_numb = uicontrol(h_fig,'Style','text','String','Number of bats','units','normalized',...
   'Position',[.01 .74 .23 .05],'fontsize',9,'fontweight','b'); 
n_bats = uicontrol(h_fig,'Style','edit','String',num2str(length(s_numbers)),'units','normalized',...
    'Position',[.01 .70 .23 .05],'fontsize',9,'fontweight','b');

%Current session
h_sess = uicontrol(h_fig,'Style','text','String','Current session','units','normalized',...
   'Position',[.26 .74 .23 .05],'fontsize',9,'fontweight','b'); 
h_session = uicontrol(h_fig,'Style','edit','String',num2str(length(session_init)),'units','normalized',...
    'Position',[.26 .70 .23 .05],'fontsize',9,'fontweight','b');

%Reward decrement
h_rew = uicontrol(h_fig,'Style','text','String','Reward decrement (s)','units','normalized',...
   'Position',[.51 .74 .23 .05],'fontsize',9,'fontweight','b'); 
h_d_rew = uicontrol(h_fig,'Style','edit','String','0.00','units','normalized',...
    'Position',[.51 .70 .23 .05],'fontsize',9,'fontweight','b');

%ExpType
h_exp = uicontrol(h_fig,'Style','text','String','Exp Type','units','normalized',...
   'Position',[.76 .74 .23 .05],'fontsize',9,'fontweight','b'); 
h_Exp = uicontrol(h_fig,'Style','popupmenu','String',{'group','solo'},'units','normalized',...
    'Position',[.76 .70 .23 .05],'fontsize',9,'fontweight','b');

%Reward types
h_rew1 = uicontrol(h_fig,'Style','text','String','Reward 1','units','normalized',...
   'Position',[.01 .64 .23 .05],'fontsize',9,'fontweight','b'); 
h_Rew1 = uicontrol(h_fig,'Style','edit','String','','units','normalized',...
    'Position',[.01 .60 .23 .05],'fontsize',9,'fontweight','b');

h_rew2 = uicontrol(h_fig,'Style','text','String','Reward 2','units','normalized',...
   'Position',[.26 .64 .23 .05],'fontsize',9,'fontweight','b'); 
h_Rew2 = uicontrol(h_fig,'Style','edit','String','','units','normalized',...
    'Position',[.26 .60 .23 .05],'fontsize',9,'fontweight','b');

h_rew3 = uicontrol(h_fig,'Style','text','String','Reward 3','units','normalized',...
   'Position',[.51 .64 .23 .05],'fontsize',9,'fontweight','b'); 
h_Rew3 = uicontrol(h_fig,'Style','edit','String','','units','normalized',...
    'Position',[.51 .60 .23 .05],'fontsize',9,'fontweight','b');

h_rew4 = uicontrol(h_fig,'Style','text','String','Reward 4','units','normalized',...
   'Position',[.76 .64 .23 .05],'fontsize',9,'fontweight','b'); 
h_Rew4 = uicontrol(h_fig,'Style','edit','String','','units','normalized',...
    'Position',[.76 .60 .23 .05],'fontsize',9,'fontweight','b');


%Feeders table
% sAlign = uistyle("HorizontalAlignment","center");
fed_table = table2cell(table((1:4)',true(4,1),true(4,1),[100;0;0;0],0.2*ones(4,1),zeros(4,1)));
h_feeders = uitable(h_fig,'Data',fed_table,'ColumnName',{'Feeder ID'; 'Enabled'; 'State'; 'Probability'; 'Duration'; 'Reward'},'units','normalized',...
    'Position',[.11 .05 .38 .50],'fontsize',9,'fontweight','b','ColumnEditable',true,'ColumnWidth',{60,50,40,65,55,60});
% addStyle(h_feeders,sAlign)

%Bat table
dim = str2double(n_bats.String);
if dim>1
    bat_table = table2cell(table(tag_names,cell2mat(bat_ids),true(dim,1),zeros(dim,1),zeros(dim,1),false(dim,1),repmat({''},dim,1)));
else
    bat_table = table2cell(table({tag_names},cell2mat(bat_ids),true(dim,1),zeros(dim,1),zeros(dim,1),false(dim,1),repmat({''},dim,1)));
end
h_bat = uitable(h_fig,'Data',bat_table,'ColumnName',{'Tag'; 'Bat ID'; 'Enabled'; 'Trials'; 'Correct'; 'Ephys'; 'Region'},'units','normalized',...
    'Position',[.50 .05 .49 .50],'fontsize',9,'fontweight','b','ColumnEditable',true,'ColumnWidth',{50,45,48,60,60,45,60});
% addStyle(h_bat,sAlign)

%feed buttons
h_feed1 = uicontrol(h_fig,'Style','pushbutton','String','Feed1','units','normalized',...
    'Position',[.01 .5 .09 .05],'Callback',@feed_button,'fontsize',10,'fontweight','b','tag','run');
h_feed2 = uicontrol(h_fig,'Style','pushbutton','String','Feed2','units','normalized',...
    'Position',[.01 .4 .09 .05],'Callback',@feed_button,'fontsize',10,'fontweight','b','tag','run');

h_feed3 = uicontrol(h_fig,'Style','pushbutton','String','Feed3','units','normalized',...
    'Position',[.01 .3 .09 .05],'Callback',@feed_button,'fontsize',10,'fontweight','b','tag','run');

h_feed4 = uicontrol(h_fig,'Style','pushbutton','String','Feed4','units','normalized',...
    'Position',[.01 .2 .09 .05],'Callback',@feed_button,'fontsize',10,'fontweight','b','tag','run');

% Run motors to change syringes
h_reload = uicontrol(h_fig,'Style','togglebutton','String','Reload','units','normalized',...
    'Position',[.01 .1 .09 .05],'Callback',@run_feeders_in_task,'fontsize',10,'fontweight','b','tag','run');


% %activate feeder toggle
% fed_override = uicontrol(h_fig,'Style','pushbutton','String','Reload','units','normalized',...
%     'Position',[.01 .1 .09 .05],'Callback',@reload_button,'fontsize',10,'fontweight','b','tag','run');


end