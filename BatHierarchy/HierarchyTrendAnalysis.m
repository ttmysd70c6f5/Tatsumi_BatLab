%% Group 1
RwrFile = dir('**\Ext_Behavior*\Feeder_Intrusion.mat');
Compute_Win_Lose(RwrFile)
p = cell(2,1);
p{1} = Compute_Win_Lose(RwrFile(1:3));
p{2} = Compute_Win_Lose(RwrFile(4:6));

%% Group 2
RwrFile = dir('**\Ext_Behavior*\Feeder_Intrusion.mat');
Compute_Win_Lose(RwrFile)
p = cell(3,1);
p{1} = Compute_Win_Lose(RwrFile(1:5));
p{2} = Compute_Win_Lose(RwrFile(6:10));
p{3} = Compute_Win_Lose(RwrFile(11:15));

%% Group 3
RwrFile = dir('**\Ext_Behavior*\Feeder_Intrusion.mat');
Compute_Win_Lose(RwrFile);
p = cell(3,1);
p{1} = Compute_Win_Lose(RwrFile(1:3));
p{2} = Compute_Win_Lose(RwrFile(4:6));
p{3} = Compute_Win_Lose(RwrFile(7:10));