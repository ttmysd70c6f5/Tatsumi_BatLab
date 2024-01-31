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
s_numbers = load('Serial_numbers_in_use.txt'); %ORDERED serial numbers!
n_bats = length(s_numbers);
client_asks = 1;
state = ones(n_bats,1);
prev_state = ones(n_bats,1);
bat_at_feeder = ones(4,2);

%[prev_state, bat_id, fed_id, client_asks];

%Before starting the communication, flush the buffer
flush(u);
while client_asks
    
    if cdp_server.BytesAvailable > 0
        flush(u);
        x = fread(cdp_server,4,'uint8');
        
        bat_id = x(2);
        prev_state(bat_id) = x(1);
        state(bat_id) = is_enabled(u,bat_id,prev_state(bat_id),s_numbers);
        fed_id = x(3);
        client_asks = x(4);
        
        tic;
        fwrite(cdp_server,uint8([state(bat_id) whos_at_feeder_v2(u,fed_id,s_numbers)]),'uint8');
        toc;
    end
     
end

fclose(cdp_server);
clear all;
