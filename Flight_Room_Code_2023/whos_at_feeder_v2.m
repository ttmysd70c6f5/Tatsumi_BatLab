function bat_n = whos_at_feeder_v2(udp_obj,feeder_n,s_n)
%Returns the id of closest bat to the input feeder
%Matrix of Feeder Coordinates
F = [2760,  900, 1650;...
    2760, -1100, 1650;...
    2760, 1450,  800;...
    2760,-1500,  800];

tags_SN = s_n;                                    %Tag Serial Numbers
n_ave = 2;                                        %Number of frames for averaging
pos = nan(length(tags_SN),3, n_ave);              %Initialize position vector

%Start reading from udp_obj (after flushing the buffer)
flush(udp_obj);
for i=1:n_ave
    counter = 0;
    while any(isnan(squeeze(pos(:,:,i))),'all') && counter < 100                    %Keep reading until all tags are acquired
        [sn,~,x,y,z] = decode_CDP_v3(udp_obj);
        if sn~=17040920                                             %Skip sync tag
            pos(find(tags_SN == sn),:,i) = [x,y,z];
        end                                                         %populate position matrix
    counter = counter + 1;
    end
end

%Evaluate closest bat
[~,bat_n] = min(vecnorm(mean(pos,3)-F(feeder_n,:),2,2));
end