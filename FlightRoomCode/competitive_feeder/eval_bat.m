function state = eval_bat(udp_obj,bat_id,prev_state,s_numbers)
% Get the new state of bat


%Evaluate if a bat is enabled to feed
sn = 0;
counter = 0;
num_buffer = 100;

if prev_state                    % If already enabled
    state = 1;
else
    %Start reading from udp_obj (after flushing the buffer)
    %Keep reading until you get the requested SN or you tried more than 100 times
    flush(udp_obj);
    while sn ~= s_numbers(bat_id) && counter < num_buffer
        [sn,~,x,~,~] = decode_CDP(udp_obj);
        counter = counter+1;
    end
    %display(counter);

    if x<700      % If crossed virtual wall
        state = 1;
    else
        state = 0;
    end
end


end
