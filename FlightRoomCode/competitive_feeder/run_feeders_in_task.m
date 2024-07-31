function run_feeders_in_task(active)
% run_feeders(active) activates the feeders specified by active. 
% The input variable active should be a scholar or a vector of feeder identities (1-4).

global h_run h_pause h_stop
global h_reload

% activate all the feeders if no feeder is specified
if isempty(who)
    active = 1 : 4;
end

speed = zeros(1,4);
speed(active) = 1; % assign a non-zero speed

run_state = logical(h_reload.Value); % 1 if button is toggled

if h_run.Value && h_pause.Value && ~h_stop.Value                            % run only during pause mode
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

    if run_state
        disp('Run all motors.');
        start(MOT1);
        start(MOT2);
        start(MOT3);
        start(MOT4);
    else
        disp('Press any key to stop all motors');
        stop(MOT1);
        stop(MOT2);
        stop(MOT3);
        stop(MOT4);
    end
end