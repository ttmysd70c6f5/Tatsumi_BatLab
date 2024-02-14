function [] = exportTrodesRecording(path_to_rec_dir, trodes_path, varargin)
%EXPORTTRODESRECORDING Summary of this function goes here
%   Detailed explanation goes here
%   Recording directory must have default name generated by trodes: 
%   (YYYYMMDD_HHMMSS.rec)
p = inputParser;
addRequired(p,'path_to_rec_dir');
addRequired(p,'path_to_trodes');

addOptional(p,'isTethered',0); % 0 for wireless recording
addOptional(p,'extractLFP',1); % 1 to extract, 0 to skip
addOptional(p,'extractAnalog',1);
addOptional(p,'extractKilosort',1);
addOptional(p,'extractTime',1);
addOptional(p,'extractDIO',1);
addOptional(p,'extractSpikes',1);
addOptional(p,'extractSpikeBand',1);
addOptional(p,'extractRaw',0);
addOptional(p,'lfpLowPass', 500);
addOptional(p,'spikesFreqBand', [600, 6000]);

parse(p,path_to_rec_dir,trodes_path,varargin{:});

isTethered = p.Results.isTethered;
extractLFP = p.Results.extractLFP;
extractAnalog = p.Results.extractAnalog;
extractKilosort = p.Results.extractKilosort;
extractTime = p.Results.extractTime;
extractDIO = p.Results.extractDIO;
extractSpikes = p.Results.extractSpikes;
extractSpikeBand = p.Results.extractSpikeBand;
extractRaw = p.Results.extractRaw;
lfpLowPass = p.Results.lfpLowPass;
spikesFreqBand = p.Results.spikesFreqBand;

fprintf("Extracting data from %s \n\n", path_to_rec_dir);

%% Select appropriate flags for extraction
possible_flags = ["-lfp", "-analogio", "-kilosort", "-dio", "-spikes", "-spikeband", "-time", "-raw"];
idx = logical([extractLFP extractAnalog extractKilosort extractDIO extractSpikes extractSpikeBand extractTime extractRaw]);
flags = possible_flags(idx);
cmd_flags = join(flags);

lfp_options = sprintf("-lfplowpass %d -lfpoutputrate 1500 -uselfpfilters 0", lfpLowPass); % LFP export parameters
spikeband_options = sprintf("-spikehighpass %d -spikelowpass %d -usespikefilters 0 -invert 0", spikesFreqBand(1), spikesFreqBand(2)); % spikes export parameters

cmd_flags = join([cmd_flags lfp_options spikeband_options]);

%% Path to Trodes is needed to call the trodes function.
addpath(genpath(trodes_path))


%% Parse recording directory
path_to_rec_file = dir(fullfile(path_to_rec_dir,'*.rec')); path_to_rec_file = fullfile(path_to_rec_file.folder,path_to_rec_file.name);
[dirpath, dirname, ext] = fileparts(path_to_rec_file);
% [dirpath, dirname, ext] = fileparts(path_to_rec_dir);
fparts = regexp(join([dirname ext],''), "(\d\d\d\d)(\d\d)(\d\d)_(\d\d\d\d\d\d)_merged.rec", "tokens");
YY = fparts{1}{1}(3:end);
MM = fparts{1}{2};
DD = fparts{1}{3};
date = sprintf("%s%s%s",YY,MM,DD);

time = fparts{1}{4};

% if(~isTethered)
%     path_to_rec_file = fullfile(path_to_rec_dir, dirname + "_merged" + ".rec");
% else
%     path_to_rec_file = fullfile(path_to_rec_dir, dirname + ".rec");
% end


original_dir = pwd;
cd(trodes_path)
if(exist(path_to_rec_file))
    formatted_path_to_rec_file = join(['"', path_to_rec_file, '"'], '');
    % output_path = join(['"', path_to_rec_dir, '"'], '');
    output_path = fullfile(path_to_rec_dir,'Processed');
    mkdir(output_path)
    
    cmd = join(['trodesexport', cmd_flags, '-rec', formatted_path_to_rec_file, '-outputdirectory', output_path, '-interp 0']); % command to run
    
    disp(cmd)
    [status, cmdout] = system(cmd); % Export lfp file
    disp(cmdout)

    
    chanmap_files = dir(join([path_to_rec_dir, '\*kilosort\*channelmap*.dat']));
    for iProbe = 1:length(chanmap_files)
        chanmap_file = chanmap_files(iProbe)
        parts = split(chanmap_file(1).name, '.');
        
        ks_struct = readTrodesExtractedDataFile(fullfile(chanmap_file(1).folder,chanmap_file(1).name));
        nchan = length(ks_struct.fields(1).data);
        chanMap = (1:1:nchan);
        chanMap0ind = chanMap-1;
        connected = logical(ones(nchan, 1));
        kcoords = ones(nchan, 1);
        xcoords = double(ks_struct.fields(1).data);
        ycoords = double(ks_struct.fields(2).data);
        fs = 30000; % sampling rate
        name = sprintf('Channel map from recording on %s_%s',date, time);
        save(fullfile(chanmap_file.folder, sprintf('%s_%s_%s.mat', parts{2}, date, time)), 'chanMap', 'chanMap0ind', 'connected','fs','kcoords','name','xcoords', 'ycoords')
        %mkdir(fullfile(path_to_rec_dir, 'kilosort_outdir'));
        %mkdir(fullfile(path_to_rec_dir, 'kilosort_workdir'));
    end
    
else
    disp("Recording file does not exist")
end
cd(original_dir);

end

