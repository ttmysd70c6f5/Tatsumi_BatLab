%% Ciholas server for Behavioral task
% This script constantly monitor the state of a tag (enabled or not)
% and calculate the closest tag to the requested feeder

%Create port for getting CDP data
u = udpport("datagram","LocalPort",7667,'EnablePortSharing',true,'Timeout',20);
configureMulticast(u,"239.255.76.67");

%Create Ciholas server to communicate with another Matlab instance
cdp_server = tcpip('0.0.0.0',30000,'NetworkRole','server');
fopen(cdp_server);

%Initialize variables
% Ordered serial numbers
% s_numbers = load('Serial_numbers_in_use.txt'); %ORDERED serial numbers!
s_table = readtable('serial_numbers.csv','Format','%d %s'); s_table = s_table(~strcmp(s_table.Bat_id,''),:);
s_numbers = double(s_table.Serial_number);

n_bats = length(s_numbers);
client_asks = 1; % communicate with CDP network while 1
% state = ones(n_bats,1);
% prev_state = ones(n_bats,1);
bat_at_feeder = ones(4,2);

%[prev_state, bat_id, fed_id, client_asks];

%Before starting the communication, flush the buffer
flush(u);

% Communicate with Foraging_Task
while client_asks
    
    if cdp_server.BytesAvailable > 0
        flush(u);
        x = fread(cdp_server,5,'uint8'); % Request from Foraging_Task: [type, bat status, bat ID, feeder ID, ask_cdp]
        request_type = x(1); % 1 = update the bat status, 2 = bat at feeder, 3 = feeder status
        client_asks = x(5); % kill the server if 0

        % 1. update bat status
        if request_type == 1
            bat_id = x(3); % Bat ID asked by client
            prev_state = x(2); % Previous bat status
            bat_state = eval_bat(u,bat_id,prev_state,s_numbers); % New bat status
            
            % Return a value to Foraging_Task
            tic;
            fwrite(cdp_server,uint8(bat_state),'uint8'); % [New state, closest bat ID]
            toc;

        % 2. identify who activated the feeder
        elseif request_type == 2
            fed_id = x(4); % Activated feeder
            bat_at_feeder = whos_at_feeder(u,fed_id,s_numbers); % distance to feeder, bat closest to the feeder
            
            % Return a value to Foraging_Task
            tic;
            fwrite(cdp_server,uint8(bat_at_feeder),'uint8'); % [New state, closest bat ID]
            toc;
        
        % 3. update the feeder status
        elseif request_type == 3
            fed_id = x(4); % Activated feeder
            prev_state = x(2); % previous feeder status
            fd_state = eval_feeder(u,fed_id,prev_state,s_numbers);
            
            % Return a value to Foraging_Task
            tic;
            fwrite(cdp_server,uint8(fd_state),'uint8'); % [New state, closest bat ID]
            toc;    
            
        end
    end
     
end

fclose(cdp_server);
clear all;
