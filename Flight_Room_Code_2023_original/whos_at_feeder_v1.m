function bat_n = whos_at_feeder_v1(udp_obj,feeder_n)
%Returns the id of closest bat to the input feeder
%Matrix of Feeder Coordinates
F = [2700,  870, 1870;...
    2700, -905, 1870;...
    2600, 1290,  900;...
    2600,-1360,  900];

tags_SN = [17106917; 17106934; 17107055; 17106969; 17107100];   %Tag Serial Numbers
n_ave = 3;                                                      %Number of frames for averaging
pos = nan(length(tags_SN),3, n_ave);                            %Initialize position vector

%Start reading from udp_obj (after flushing the buffer)
flush(udp_obj);
for i=1:n_ave
    while any(isnan(squeeze(pos(:,:,i))),'all')                                     %Keep reading until all tags are acquired
        [sn,~,x,y,z] = decode_CDP_v3(udp_obj);
        if sn~=17106963                                             %Skip sync tag
            pos(find(tags_SN == sn),:,i) = [x,y,z];
        end                                                         %populate position matrix
    end
end

%Evaluate closest bat
[~,bat_n] = min(vecnorm(mean(pos,3)-F(feeder_n,:),2,2));
end