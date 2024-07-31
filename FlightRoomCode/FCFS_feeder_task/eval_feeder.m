function state = eval_feeder(udp_obj,feeder_n,prev_state,s_numbers)
% Returns the new state of the feeder
% udp_obj: udp object
% feeder_n: feeder id
% prev_state: previous state of the feeder
% s_n: list of serial numbers

d_th = 300; % Threshold 
num_buffer = 100; % Buffer size

%Matrix of Feeder Coordinates
F = [2760,  900, 1650;...
    2760, -1100, 1650;...
    2760, 1450,  800;...
    2760,-1500,  800];

n_ave = 2;                                          %Number of frames for averaging
pos = nan(length(s_numbers),3, n_ave);              %Initialize position vector [sn,dim,frame]

if prev_state == 1                          % If already enabled
    state = 1;
else
    %Start reading from udp_obj (after flushing the buffer)
    flush(udp_obj);
    for i=1:n_ave
        counter = 0;
        while any(isnan(squeeze(pos(:,:,i))),'all') && counter < num_buffer                    %Keep reading until all tags are acquired
            [sn,~,x,y,z] = decode_CDP(udp_obj); % [serial number,network time,x,y,z] 
            if sn~=17040920                                             %Skip sync tag
                pos(find(s_numbers == sn),:,i) = [x,y,z];
            end                                                         %populate position matrix
        counter = counter + 1;
        end
    end
    
    %Evaluate the distance to feeder for all bats
    d_feeder = vecnorm(mean(pos,3)-F(feeder_n,:),2,2); % distances to feeder: [sn,1]
    if all(d_feeder(~isnan(d_feeder)) > d_th) % if the feeder is emptied
        state = 1;
    else
        state = 0;
    end
end