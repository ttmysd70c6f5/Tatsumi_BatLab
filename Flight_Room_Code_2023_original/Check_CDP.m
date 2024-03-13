%% Check the serial numbers of the availbale tags

%If you are not getting packets, try to restart this computer WHILE the
%Ciholas is broadcasting data!!!


%Create port for getting CDP data
u = udpport("datagram","LocalPort",7667,'EnablePortSharing',true);
configureMulticast(u,"239.255.76.67");

for i = 1:1

s = 0;
sn = [];
counter = 0;

flush(u);
tic
while s<8 && counter < 500
    [test_SN,~,~,~,~] = decode_CDP_v3(u);
    sn = [sn test_SN];
    s = length(unique(sn));
    [C,ia,ic] = unique(sn);
    counter = counter +1;
end
toc

attempt(i) = counter;
end

histogram(attempt);
disp(unique(sn));
close all;