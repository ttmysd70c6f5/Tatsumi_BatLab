function state = is_enabled(udp_obj,bat_id,pre_state,s_n)
%Evaluate if a bat is enabled to feed
sn = 0;
counter = 0;
tags_SN = s_n;

%Start reading from udp_obj (after flushing the buffer)
%Keep reading until you get the requested SN or you tried more than 100 times
flush(udp_obj);
while sn ~= tags_SN(bat_id) && counter < 100
    [sn,~,x,~,~] = decode_CDP_v3(udp_obj);
    counter = counter+1;
end
%display(counter);
if pre_state || x<700      %If already enabled or crossed virtual wall
    state = 1;
else
    state = 0;
end
end

