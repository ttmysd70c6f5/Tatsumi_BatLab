% Load data
% Reward signals
RewardFile = dir('Ext_Ephys*/*c3d*.mat');
load(fullfile(RewardFile(1).folder,RewardFile(1).name),'AnalogFrameRate','AnalogSignals') % reward signals from the feeder (1st col)
% Behavior
BehFile = dir('Ext_Behavior*\Extracted_Behavior*.mat');
load(fullfile(BehFile.folder,BehFile.name))

% Behavior data

t_rwrd = find(diff(AnalogSignals(:,1))>2);
t_rwrd = t_rwrd - find(diff(AnalogSignals(:,2))>2,1,'first');
t_rwrd = t_rwrd / AnalogFrameRate; % convert to sec
t_rwrd = round(t_rwrd*2)/2;

rwrd = zeros(T,1);
for tt = 1:length(t_rwrd)
    rwrd(t == t_rwrd(tt)) = 1;
end

figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0.05 0.9 0.85]);
tiledlayout(4, 1, 'TileSpacing', 'tight');

for i=1:4
    ax(i) = nexttile;
    hold on
    plot(t,bflying(:,i)/2 ,'b')
    plot(t,rwrd,'r')
    hold off
    
end
linkaxes(ax)