% directory path
ephys_dir{1} = 'Z:\users\Tatsumi\data\MultiReward\multibat\240208\ephys\20240208_132817.rec\14633'; % ephys directory containing the merged data for the data 1
ephys_dir{2} = 'Z:\users\Tatsumi\data\MultiReward\singlebat\14633\240208\ephys\20240208_082558.rec\14633';
outdir = 'Z:\users\Tatsumi\data\MultiReward\multi_single_concat\240208\ephys\concatenated.rec\14633\Processed'; % path to the directory to output the concatenated data
trodes_path = 'C:\Users\Tatsumi\Documents\Trodes\Trodes_2-5-0_Windows64'; % trodes 2.5.0

if(~exist(outdir));  mkdir(outdir);  end

% setting parameters
isTethered = false; % true if tetherd recording

% Create the command to extract Trodes
for i = 1:length(ephys_dir)
    if(~isTethered)
        recDir = dir(fullfile(ephys_dir{i},'*merged.rec'));
    else
        recDir = dir(fullfile(ephys_dir{i},'*.rec'));
    end

    path_to_rec_file{i} = fullfile(recDir(1).folder,recDir(1).name);

    [dirpath, dirname, ext] = fileparts(path_to_rec_file{i});
    fparts = regexp(join([dirname ext],''), "(\d\d\d\d)(\d\d)(\d\d)_(\d\d\d\d\d\d).*.rec", "tokens");

    YY = fparts{1}{1}(3:end);
    MM = fparts{1}{2};
    DD = fparts{1}{3};
    time = fparts{1}{4};
    date = sprintf("%s%s%s",YY,MM,DD);
end

original_dir = pwd;
cd(trodes_path)

% export_flags = ["-lfp", "-analogio", "-kilosort", "-dio", "-spikes", "-spikeband", "-time", "-raw"];
export_flags = ["-lfp", "-analogio", "-kilosort", "-dio", "-spikeband", "-time"];

formatted_flags = strjoin(cellfun(@(x) x, export_flags, 'UniformOutput', false));
formatted_path_to_rec_file = strjoin(cellfun(@(x) join(['-rec', ' "', x, '"']), path_to_rec_file, 'UniformOutput', false));
out_path = ['"', outdir, '"'];

lfp_flag = '-lfplowpass 500 -lfpoutputrate 1500 -uselfpfilters 0'; % LFP export parameters
spikes_flag = '-spikehighpass 600 -spikelowpass 6000 -usespikefilters 0'; % spikes export parameters

cmd = join(['trodesexport', export_flags, formatted_path_to_rec_file, '-outputdirectory', out_path, '-interp 0', spikes_flag, lfp_flag]);

%%
% Run export function
[status, cmdout] = system(cmd); % Export lfp file
disp(cmdout)

% Create channel map for kilosort
chanmap_files = dir(join([out_path, '\*kilosort\*channelmap*.dat']));
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
    
cd(original_dir);


%% Parameters
% dates = ["231007", "231008"]
% batID = "32622"
% isTethered = false;
% 
% %outPath = fullfile(dataPath, batID, strcat(experiments(1), "_", experiments(2)), 'raw', date, 'ephys');
% outPath = fullfile(dataPath, batID, experiments(1), 'raw', strcat(dates(1), "_", dates(2)), 'ephys');
% 

%% Experiment specific parameters (should be the same throughout experiment)
% dataPath = 'Z:\users\KevinQi\datasets\GridBat';
% 
% trodesPath = 'Z:\users\KevinQi\Trodes_2-4-1';
% addpath(genpath('C:\Users\batlab\Desktop\pixel-bat'));
% experiments = ["flight_room"];%, "big_cage"];%, 'small_cage'];
% idx = 1;
% 
% for i = 1:length(dates)
%     for iExp = 1:length(experiments)
%         rawPaths{idx} = fullfile(dataPath, batID, experiments(iExp), 'raw', date);
%         processedPaths{idx} = fullfile(dataPath, batID, experiments(iExp), 'processed', date);
%         ephysPaths{idx} = fullfile(dataPath, batID, experiments(iExp), 'raw', date, 'ephys');
%         idx = idx + 1;
%     end
% end
% 
% if(~exist(outPath))
%     mkdir(outPath)
% end

%% Parse recording directory
% for i = 1:length(ephysPaths)
%     recs = dir(fullfile(ephysPaths{i}, "*.rec")); % list of recording directory
%     path_to_rec_dir{i} = fullfile(recs(1).folder, recs(1).name);
% 
%     [dirpath, dirname, ext] = fileparts(path_to_rec_dir{i});
%     fparts = regexp([dirname ext], "(\d\d\d\d)(\d\d)(\d\d)_(\d\d\d\d\d\d)*.rec", "tokens");
%     YY = fparts{1}{1}(3:end);
%     MM = fparts{1}{2};
%     DD = fparts{1}{3};
%     date = sprintf("%s%s%s",YY,MM,DD);
% 
%     time = fparts{1}{4};
% 
%     % Export LFP, Spikes, and Kilosort data
%     lfp_flag = '-lfplowpass 500 -lfpoutputrate 1500 -uselfpfilters 0'; % LFP export parameters
%     spikes_flag = '-spikehighpass 600 -spikelowpass 6000 -usespikefilters 0'; % spikes export parameters
%     if(~isTethered)
%         path_to_rec_file(i) = fullfile(path_to_rec_dir{i}, dirname + "_merged" + ".rec");
%     else
%         path_to_rec_file(i) = fullfile(path_to_rec_dir{i}, dirname + ".rec");
%     end
% end
% 
% original_dir = pwd;
% cd(trodesPath)
% formatted_path_to_rec_file_1 = join(['"', path_to_rec_file(1), '"'], '');
% formatted_path_to_rec_file_2 = join(['"', path_to_rec_file(2), '"'], '');
% out_path = join(['"', outPath, '"'], '');
% 
% cmd = join(['trodesexport -lfp -spikes -spikeband -time -dio -analogio -kilosort -rec', formatted_path_to_rec_file_1, '-rec', formatted_path_to_rec_file_2, '-outputdirectory', out_path, '-interp 0', spikes_flag lfp_flag]); % command to run
% [status, cmdout] = system(cmd); % Export lfp file
% disp(cmdout)
% 
% 
% chanmap_files = dir(join([outPath, '\*kilosort\*channelmap*.dat']));
% for iProbe = 1:length(chanmap_files)
%     chanmap_file = chanmap_files(iProbe)
%     parts = split(chanmap_file(1).name, '.');
% 
%     ks_struct = readTrodesExtractedDataFile(fullfile(chanmap_file(1).folder,chanmap_file(1).name));
%     nchan = length(ks_struct.fields(1).data);
%     chanMap = (1:1:nchan);
%     chanMap0ind = chanMap-1;
%     connected = logical(ones(nchan, 1));
%     kcoords = ones(nchan, 1);
%     xcoords = double(ks_struct.fields(1).data);
%     ycoords = double(ks_struct.fields(2).data);
%     fs = 30000; % sampling rate
%     name = sprintf('Channel map from recording on %s_%s',date, time);
%     save(fullfile(chanmap_file.folder, sprintf('%s_%s_%s.mat', parts{2}, date, time)), 'chanMap', 'chanMap0ind', 'connected','fs','kcoords','name','xcoords', 'ycoords')
%     %mkdir(fullfile(path_to_rec_dir, 'kilosort_outdir'));
%     %mkdir(fullfile(path_to_rec_dir, 'kilosort_workdir'));
% end
% 
% cd(original_dir);


