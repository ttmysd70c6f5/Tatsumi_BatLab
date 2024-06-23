function run_feeders(active)
% run_feeders(active) activates the feeders specified by active. 
% The input variable active should be a scholar or a vector of feeder identities (1-4).

% activate all the feeders if no feeder is specified
if isempty(who)
    active = 1 : 4;
end

speed = zeros(1,4);
speed(active) = 1; % assign a non-zero speed

% Connect to arduino
Port ='com3';
ID = 'uno';
ard=arduino(Port,ID,'Libraries','Adafruit\MotorShieldV2'); % Connect to Arduino hardware
shield=addon(ard,'Adafruit\MotorShieldV2'); % create a connection to an Adafruit MotorShield
% creates a connection to the DC motor
MOT1=dcmotor(shield, 1);
MOT2=dcmotor(shield, 2);
MOT3=dcmotor(shield, 3);
MOT4=dcmotor(shield, 4);
% set a speed to the DC motor
MOT1.Speed = speed(1);
MOT2.Speed = speed(2);
MOT3.Speed = speed(3);
MOT4.Speed = speed(4);

disp('press any key to run all motors');
pause; 
start(MOT1);
start(MOT2);
start(MOT3);
start(MOT4);
disp('press any key to stop all motors');
pause;
stop(MOT1);
stop(MOT2);
stop(MOT3);
stop(MOT4);
