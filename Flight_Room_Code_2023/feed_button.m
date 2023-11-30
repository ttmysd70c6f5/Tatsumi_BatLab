function feed_button(varargin)

global rewardspeed 
global h_feed1 h_feed2 h_feed3 h_feed4 
global feed1 feed2 feed3 feed4

%set feeder speeds
feed_speed = rewardspeed;
feed1.Speed = feed_speed; feed2.Speed = feed_speed; feed3.Speed = feed_speed; feed4.Speed = feed_speed;
rew_dur = 0.2;      %default value

%determine which feeder was activated
activated_feeder = find([h_feed1.Value h_feed2.Value h_feed3.Value h_feed4.Value]);
feeder_name = ['feed' num2str(activated_feeder)];

start(eval(feeder_name));
pause(rew_dur);
stop(eval(feeder_name));

%turn off all toggle buttons
h_feed1.Value = 0; h_feed2.Value = 0; h_feed3.Value = 0; h_feed4.Value = 0;