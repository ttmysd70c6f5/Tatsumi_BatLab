function reload_button(varargin)

global fed_override
global rewardspeed 
global feed1 feed2 feed3 feed4

%set feeder speeds
feed_speed = rewardspeed;
feed1.Speed = feed_speed; feed2.Speed = feed_speed; feed3.Speed = feed_speed; feed4.Speed = feed_speed;

start(eval('feed1'));
disp('Loading Feeder...');
pause(15);                  %% Make sure this is shorter than the UDP timeout (or modify that)
stop(eval('feed1'));
disp('Done Loading Feeder...');

fed_override.Value = 0;
end