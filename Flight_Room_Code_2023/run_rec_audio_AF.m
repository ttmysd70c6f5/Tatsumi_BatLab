function run_rec_audio_AF(varargin)

%records audio and sync pulses
global h_directory_aud h_run_aud h_stop_aud
global IDsound fs rec_dur input_channels inter_file_time
global batName dateSesh

timeSesh = datestr(now,'hhMMss');
homeDir_full = [h_directory_aud.String '\' timeSesh];
if ~isdir(homeDir_full)
    mkdir(homeDir_full)
end

%prepare soundmex
if 1 ~= soundmexpro('init','driver',IDsound,'samplerate',fs,'input',input_channels)
    error(['error calling ''init''' error_loc(dbstack)]);
end

soundmexpro('recbufsize', 'value', rec_dur*fs);
soundmexpro('recpause','value', ones(1,length(input_channels)),'channel',input_channels);

soundmexpro('start','length', 0)

[ret,fs,bufsiz]=soundmexpro('getproperties');
filecounter = 0;
waitless = 0;
pause(10) %small pause to make sure the first file is right length
while h_stop_aud.Value < 1 
    tic
    pause(inter_file_time-waitless)
    filecounter = filecounter + 1;
    %get audio (should include ttl)
    %[success, recbuf, pos] = soundmexpro('recgetdata', 'channel',1:7);     %This was the previous version
    [success, recbuf, pos] = soundmexpro('recgetdata', 'channel',[0:3,7]);  % AF 2023, to record only relevant mics (1 to 4) + TTL
    clocktime = clock;
    audiofile = [homeDir_full filesep batName '_' dateSesh '_audio_trial_' num2str(filecounter)]
    audtic = tic;
    %save audio
    save(audiofile, 'recbuf','fs','bufsiz','audiofile','clocktime')
    waitless = toc(audtic)
    toc
end

soundmexpro('stop');
soundmexpro('exit');

%reset the values of the run and stop toggle buttons to zero
h_run_aud.Value = 0;
h_stop_aud.Value = 0;








