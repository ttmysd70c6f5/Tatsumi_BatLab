function [digitalData] = loadTrodesDigital(path_to_recording_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% [dirpath, dirname, ext] = fileparts(path_to_recording_dir);
% assert(strcmp(ext, ".rec"), "Trodes recording directory must end in .rec");

% mergedDigital_dirname = dirname + "_merged.DIO";

mergedDigital_dir = dir(fullfile(path_to_recording_dir,'*.DIO','*dio_Controller_Din*.dat'));

%% Load data from 6 digital input channels
for i = 1:length(mergedDigital_dir)
    % mergedDigital_filename = dirname + sprintf("_merged.dio_Controller_Din%d.dat", i);
    mergedDigital_path = fullfile(mergedDigital_dir(i).folder,mergedDigital_dir(i).name);
    
    % data = readTrodesExtractedDataFile(fullfile(path_to_recording_dir, ...
                                                % mergedDigital_dirname, ...
                                                % mergedDigital_filename));
    data = readTrodesExtractedDataFile(mergedDigital_path);
    digitalState = data.fields(:,2).data;
    digitalTime = data.fields(:,1).data;

    digitalData{i}.state = digitalState;
    digitalData{i}.stateSampleTimes = digitalTime;
    digitalData{i}.raw = data;
    digitalData{i}.clockrate = data.clockrate;
    digitalData{i}.first_timestamp = data.first_timestamp;

    digitalData{i}.ttl_timestamp_usec = double(digitalTime((digitalState == 1)))*1e6/data.clockrate;
    digitalData{i}.first_timestamp_usec = double(data.first_timestamp)*1e6/data.clockrate;
end


end