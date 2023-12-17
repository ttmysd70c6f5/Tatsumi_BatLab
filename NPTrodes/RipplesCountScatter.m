function RipplesCountScatter(lfp_dir)
%% Experiment info
gain = 50; % LFP gain
fs = 1500; % LFP sampling rate
nChan = 384;

%% Load data
lfp_file = dir(fullfile(lfp_dir,'*LFP_nt*.dat')); % list of LFP .dat files

ch_depth = zeros(nChan,1);
ripple_cnt = zeros(nChan,1);

for ch = 1:nChan
    fprintf('Extracting the ripple for electrode %d\n',ch)
    lfp_data = readTrodesExtractedDataFile(fullfile(lfp_file(ch).folder, lfp_file(ch).name)); % Load LFP data
    lfp = double(lfp_data.fields.data);
    str = extract(lfp_file(ch).name,digitsPattern);
    ch_depth(ch) = 10*(str2double(str{3})-1000);
    
    % Subtract median
    load(fullfile(lfp_file(1).folder,'lfpmedian.mat'),'lfpmed1')
    lfp = lfp - lfpmed1;
    lfp = lfp - median(lfp);
    
    lfp = lowpass(lfp,400,fs);
    
    % lfp = lfp / gain; % convert to µV
    
    % Extract ripples
    Ripples = ExtractRipples(lfp,false);
    ripple_cnt(ch) = length(Ripples);
end

%% FIGURE
% scatter(ch_depth,log10(ripple_cnt),'k.')
scatter(ch_depth,ripple_cnt,'k.')
xlabel('Depth from tip (um)')
ylabel('Count of detected ripples')