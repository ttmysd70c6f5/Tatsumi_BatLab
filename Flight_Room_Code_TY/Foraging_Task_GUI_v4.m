function Foraging_Task_GUI_v4
clc;
clear all;
clear global;

%Load configuration file
load('FR_config.mat');

%Declare global variables
global h_directory h_run h_stop n_bats h_bat h_feeders
global lightbar1 lightbar2 lightbar3 lightbar4
global feed1 feed2 feed3 feed4 ard trg
global rewardspeed h_d_rew
global h_feed1 h_feed2 h_feed3 h_feed4 
global s_numbers;
global fed_override;

%Insert here the serial numbers of the tags (ORDERED)
s_numbers = load('Serial_numbers_in_use.txt');
bat_names = ['Bt1'; 'Bt2'; 'Bt3'; 'Bt4'; 'Bt5'; 'Bt6'; 'Bt7'];
bat_names = bat_names(1:length(s_numbers),:);
rewardspeed = 1; 

%Initialize and configure arduino
ard=arduino('com3','uno','Libraries','Adafruit\MotorShieldV2');
shield=addon(ard,'Adafruit\MotorShieldV2');
trg = 'D7';         configurePin(ard,trg,'pullup');     configurePin(ard,trg,'DigitalOutput');      writeDigitalPin(ard,trg,0);
lightbar1 = 'D3';   configurePin(ard,lightbar1,'pullup');
lightbar2 = 'D2';   configurePin(ard,lightbar2,'pullup');
lightbar3 = 'D4';   configurePin(ard,lightbar3,'pullup');
lightbar4 = 'D5';   configurePin(ard,lightbar4,'pullup');
MOT1 = 2;           feed1=dcmotor(shield,MOT1);     feed1.Speed = rewardspeed;
MOT2 = 1;           feed2=dcmotor(shield,MOT2);     feed2.Speed = rewardspeed; 
MOT3 = 3;           feed3=dcmotor(shield,MOT3);     feed3.Speed = rewardspeed;
MOT4 = 4;           feed4=dcmotor(shield,MOT4);     feed4.Speed = rewardspeed;
disp('--ARDUINO correctly initialized--');

%Set seed for random number generation
rand('state',sum(100*clock)); 

%% Create GUI
h_fig = figure('name','4 Feeders task','Position',[65 143 700 500],'ToolBar','none','DockControls','on','MenuBar','none');  

%default directory
h_dir = uicontrol(h_fig,'Style','text','String','Directory and File Name','units','normalized',...
   'Position',[.01 .84 .3 .05],'fontsize',9,'fontweight','b'); 
h_directory = uicontrol(h_fig,'Style','edit','String',[cd '\Data\' datestr(now,'mmddyy_HH_MM')],'units','normalized',...
    'Position',[.01 .80 .96 .05],'fontsize',9,'fontweight','b');

%Number of bats
h_numb = uicontrol(h_fig,'Style','text','String','Number of bats','units','normalized',...
   'Position',[.01 .74 .3 .05],'fontsize',9,'fontweight','b'); 
n_bats = uicontrol(h_fig,'Style','edit','String',num2str(length(s_numbers)),'units','normalized',...
    'Position',[.01 .70 .3 .05],'fontsize',9,'fontweight','b');

%Reward decrement
h_rew = uicontrol(h_fig,'Style','text','String','Reward decrement (s)','units','normalized',...
   'Position',[.4 .74 .3 .05],'fontsize',9,'fontweight','b'); 
h_d_rew = uicontrol(h_fig,'Style','edit','String','0.00','units','normalized',...
    'Position',[.4 .70 .3 .05],'fontsize',9,'fontweight','b');

%hitting run starts the task
h_run = uicontrol(h_fig,'Style','togglebutton','String','Run','units','normalized',...
    'Position',[.01 .9 .3 .1],'Callback',@Foraging_Task_v3,'fontsize',10,'fontweight','b','tag','run');

%stop task
h_stop = uicontrol(h_fig,'Style','togglebutton','String','Stop','units','normalized',...
    'Position',[.32 .9 .3 .1],'fontsize',10,'fontweight','b','tag','run');

%Feeders table
fed_table = table2cell(table((1:4)',[100;0;0;0],0.2*ones(4,1)));
h_feeders = uitable(h_fig,'Data',fed_table,'ColumnName',{'Feeder ID'; 'Probability'; 'Duration'},'units','normalized',...
    'Position',[.11 .15 .38 .50],'fontsize',9,'fontweight','b','ColumnEditable',true);

%Bat table
dim = str2double(n_bats.String);
if dim>1
    bat_table = table2cell(table(bat_names,true(dim,1),zeros(dim,1),zeros(dim,1)));
else
    bat_table = table2cell(table({bat_names},true(dim,1),zeros(dim,1),zeros(dim,1)));
end
h_bat = uitable(h_fig,'Data',bat_table,'ColumnName',{'Bat ID'; 'Enabled'; 'Trials'; 'Correct';},'units','normalized',...
    'Position',[.50 .15 .50 .50],'fontsize',9,'fontweight','b','ColumnEditable',true);

%feed buttons
h_feed1 = uicontrol(h_fig,'Style','pushbutton','String','Feed1','units','normalized',...
    'Position',[.01 .5 .09 .05],'Callback',@feed_button,'fontsize',10,'fontweight','b','tag','run');
h_feed2 = uicontrol(h_fig,'Style','pushbutton','String','Feed2','units','normalized',...
    'Position',[.01 .4 .09 .05],'Callback',@feed_button,'fontsize',10,'fontweight','b','tag','run');

h_feed3 = uicontrol(h_fig,'Style','pushbutton','String','Feed3','units','normalized',...
    'Position',[.01 .3 .09 .05],'Callback',@feed_button,'fontsize',10,'fontweight','b','tag','run');

h_feed4 = uicontrol(h_fig,'Style','pushbutton','String','Feed4','units','normalized',...
    'Position',[.01 .2 .09 .05],'Callback',@feed_button,'fontsize',10,'fontweight','b','tag','run');

%activate feeder toggle
fed_override = uicontrol(h_fig,'Style','pushbutton','String','Reload','units','normalized',...
    'Position',[.01 .1 .09 .05],'Callback',@reload_button,'fontsize',10,'fontweight','b','tag','run');


end