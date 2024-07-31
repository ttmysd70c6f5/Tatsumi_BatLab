function test_beam_breaks()
% RUN THIS TEST WITH LOW LIGHT LEVELS!!
Port ='com3';
ID = 'uno';
ard=arduino(Port,ID,'Libraries','Adafruit\MotorShieldV2');
shield=addon(ard,'Adafruit\MotorShieldV2');
feed1=dcmotor(shield, 2);    feed1.Speed = 1; lightbar1 = 'D3';   configurePin(ard,lightbar1,'pullup');
feed2=dcmotor(shield, 1);    feed2.Speed = 1; lightbar2 = 'D2';   configurePin(ard,lightbar2,'pullup');
feed3=dcmotor(shield, 3);    feed3.Speed = 1; lightbar3 = 'D4';   configurePin(ard,lightbar3,'pullup');
feed4=dcmotor(shield, 4);    feed4.Speed = 1; lightbar4 = 'D5';   configurePin(ard,lightbar4,'pullup');

while 1
    
    %Reset state of beams
    lb = ones(4,1);
    
    %Check for beam breaks
    while all(lb)
        lb(1) = readDigitalPin(ard,lightbar1);
        lb(2) = readDigitalPin(ard,lightbar2);
        lb(3) = readDigitalPin(ard,lightbar3);
        lb(4) = readDigitalPin(ard,lightbar4);
        disp('...');
    end
    
    %If a beam break occurs
    if ~all(lb)
        interrupted_beams = find(~lb);
        for bb = 1:numel(interrupted_beams)
            active_feeder = interrupted_beams(bb);
            feeder_name = ['feed' num2str(active_feeder)];
            start(eval(feeder_name));   pause(0.1);  stop(eval(feeder_name));
            disp(feeder_name);
        end
    end
end
end

