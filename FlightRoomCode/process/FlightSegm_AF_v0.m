function [f_vector,f_num,f_str,f_stp] = FlightSegm_AF_v0(v_abs,v_th,Fs)
%%Detect flight starts and stops by using risetime and falltime
% ---- OUTPUT ----
% f_vector = column vector of 0s (not flying) and 1s (flying)
% f_num = number of complete flights
% f_str = samples corresponding to take-off
% f_stp = samples corresponding to landing

% Initialize vectors and detect flights
f_vector = zeros(size(v_abs));
allsums = normalize(movsum(v_abs > v_th,Fs,1),'range');

% Detect take-off and landing. More than 40% lower level reccommended
[~,rLT,~,~,~] = risetime(allsums,'PercentReferenceLevels',[50 95]);    f_str = round(rLT);
[~,fLT,~,~,~] = falltime(allsums,'PercentReferenceLevels',[50 95]);    f_stp = round(fLT);

% Ensure same number of starts/stops
if ~isempty(f_str)
    f_str = f_str(f_str<f_stp(end,1));
    f_stp = f_stp(f_stp>f_str(1,1));
    f_num = size(f_stp,1);
else
    f_str = [];
    f_stp = [];
    f_num = 0;
end

% Fill-up the f_vector
for f = 1:f_num
    f_vector(f_str(f,1):f_stp(f,1),1) = 1;
end

end
