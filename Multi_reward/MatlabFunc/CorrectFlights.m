function Flights = CorrectFlights(Flights,t_start,t_end)
n_bats = length(Flights);

for flying_bat = 1:n_bats
    if ~isempty(Flights(flying_bat).t_sec)
        t_sec = Flights(flying_bat).t_sec;
        n_flight = Flights(flying_bat).n_flight;

        % Exclude the flight at very biggining or end
        is_include = t_sec > t_start & t_sec < t_end;
        

        Flights(flying_bat).is_include = is_include;
    end
end

% function Flights_corrected = CorrectFlights(BhvData,Flights,prms)
% 
% t_margin = prms.t_margin;
% t_ephys_start = prms.t_ephys_start;
% t_ephys_end = prms.t_ephys_end;
% n_bats = length(Flights);
% t = BhvData.t;
% 
% for flying_bat = 1:n_bats
%     if ~isempty(Flights(flying_bat).t_sec)
%         t_sec = Flights(flying_bat).t_sec;
% 
%         % Exclude the flight at very biggining or end
%         is_include_edge = t_sec > t(1)-t_margin(1) & t_sec < t(end)-t_margin(2);
% 
%         % Exclude the flight without ephys recording
%         is_include_ephys = t_sec > t_ephys_start & t_sec < t_ephys_end;
% 
%         idx = is_include_edge & is_include_ephys; % index to remmain
% 
%         Flights_corrected(flying_bat).t_idx = Flights(flying_bat).t_idx(idx);
%         Flights_corrected(flying_bat).t_sec = Flights(flying_bat).t_sec(idx);
%         Flights_corrected(flying_bat).a_abs_observe = Flights(flying_bat).a_abs_observe(idx);
%         Flights_corrected(flying_bat).clu_idx = Flights(flying_bat).clu_idx(idx);
%     end
% end
% 
% 

