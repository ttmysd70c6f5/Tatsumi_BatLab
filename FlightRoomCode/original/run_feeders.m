
function run_feeders(active)

%enter feeders to activate

if isempty(who)
    active = 1 : 4;
end

speed = zeros(1,4);
speed(active) = 1;

Port ='com3';
ID = 'uno';
ard=arduino(Port,ID,'Libraries','Adafruit\MotorShieldV2');
shield=addon(ard,'Adafruit\MotorShieldV2');
MOT1=dcmotor(shield, 2);
MOT2=dcmotor(shield, 1);
MOT3=dcmotor(shield, 3);
MOT4=dcmotor(shield, 4);

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
