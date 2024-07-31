function run_feeders_in_task(varargin)
% run_feeders(active) activates the feeders specified by active. 
% The input variable active should be a scholar or a vector of feeder identities (1-4).

global h_run h_pause h_stop
global h_reload
global feed1 feed2 feed3 feed4

run_state = logical(h_reload.Value); % 1 if button is toggled

if h_run.Value && h_pause.Value && ~h_stop.Value                            % run only during pause mode
    if run_state        % Run feeder
        start(feed1)
        start(feed2)
        start(feed3)
        start(feed4)
    else
        stop(feed1)
        stop(feed2)
        stop(feed3)
        stop(feed4)
    end
end