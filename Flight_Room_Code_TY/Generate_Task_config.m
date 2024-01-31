% Script to generate configuration file for Flight Room behavior

num_bats = 5;
tags_SN = [17106917; 17106934; 17107055; 17106969; 17107100];
bat_nms = ['Dai'; 'Den'; 'Dia'; 'Dor'; 'Dum'];
F = [2700,  870, 1870;...
     2700, -905, 1870;...
     2600, 1290,  900;...
     2600,-1360,  900];
                
save('FR_config.mat');

