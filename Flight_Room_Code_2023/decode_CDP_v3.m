function [sn,nt,x,y,z] = decode_CDP_v3(udp_obj,~)
%A CDP Packet is made up of a CDP Packet Header followed a list of CDP Data Items.
type = 0;

%Only get position V3 data
while type ~= 309
    
    packet = read(udp_obj,1);
    data = uint8(packet.Data);
    
    %Cut CDP Packet Header
    data = data(21:end);
    
    %After decoding the Heather, look for CDP Data Items: Type, Size and Actual Data
    type = typecast(data(1:2),'uint16');               data = data(3:end);
    size = typecast(data(1:2),'uint16');               data = data(3:end);
    di_data = data(1:size);
    
    if type == 309
    sn_p = typecast(di_data(1:4),'uint32');                        di_data = di_data(5:end);
    nt_p = double(typecast(di_data(1:8),'int64'))*15.65e-12;       di_data = di_data(9:end);
    x_p =  typecast(di_data(1:4),'int32');                         di_data = di_data(5:end);
    y_p =  typecast(di_data(1:4),'int32');                         di_data = di_data(5:end);
    z_p =  typecast(di_data(1:4),'int32');                         di_data = di_data(5:end);
    end
end
sn = sn_p;   nt = nt_p;    x = x_p;    y = y_p;   z= z_p;
end