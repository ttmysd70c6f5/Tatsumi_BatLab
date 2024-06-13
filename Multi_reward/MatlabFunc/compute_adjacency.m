function varargout = compute_adjacency(BhvData,varargin)

p = inputParser;
addRequired(p, 'BhvData');
addOptional(p, 't_range', []); % all = all landings, fd = only feeder landings
addOptional(p, 'shuffle', false);
addOptional(p, 'type', 'spatial');

parse(p, BhvData, varargin{:});
t_range = p.Results.t_range;
if_shuffle = p.Results.shuffle;
adj_type = p.Results.type;

% Load features
blanding = 1-BhvData.bflying;
t = BhvData.t;
T = BhvData.T;
r = BhvData.r;
n_bats = BhvData.n_tags;

if ~isempty(t_range)
    idx = t >= t_range(1) & t <= t_range(2);
    t = t(idx);
    r = r(idx,:,:);
    blanding = blanding(idx,:);
end

if if_shuffle
    shft = randi(T,[n_bats,n_bats]);
end

% Threshold for adjacency evaluation
th_adj = 0.5;
coeff_adj = th_adj^2/(log(2)); % scaling coefficient

% Identify when all bats are resting
% blanding_all = all(blanding,2);

% Create adjacency matrix
switch adj_type
    case 'spatial'
        adj_trace = cell(n_bats,n_bats);
        t_adj_trace = cell(n_bats,n_bats);
        len_adj_trace = cell(n_bats,n_bats);
    case 'time'
        adj_mat = NaN(n_bats,n_bats);
end

% Specify a pair of bat
for bb = 1:n_bats
    for bbb = 1:n_bats
        if bb ~= bbb
            if if_shuffle
                blanding_all = blanding(:,bb) .* circshift(blanding(:,bbb),shft(bb,bbb)); % 1 when both pairwise bats are resting
                r = circshift(r,shft(bb,bbb),1);
            else
                blanding_all = blanding(:,bb) .* blanding(:,bbb); % 1 when both pairwise bats are resting
            end
            
            % Identify the time when new states happen
            t_start_idx = find(diff([blanding_all(1); blanding_all])==1); % timestamp when the network changes
            t_end_idx = find(diff([blanding_all(1); blanding_all])==-1); % timestamp when the network changes
            if blanding_all(1) == 0 && blanding_all(end) == 0
            elseif blanding_all(1) == 0 && blanding_all(end) == 1
                t_start_idx = t_start_idx(1:end-1);
            elseif blanding_all(1) == 1 && blanding_all(end) == 0
                t_end_idx = t_end_idx(2:end);
            elseif blanding_all(1) == 1 && blanding_all(end) == 1
                t_start_idx = t_start_idx(1:end-1);
                t_end_idx = t_end_idx(2:end);
            end

            t_start_sec = t(t_start_idx); % time when the network changes
            t_end_sec = t(t_end_idx);
            
            t_len_sec = t_end_sec-t_start_sec; % length of each state
            t_sum_sec = sum(t_len_sec);

            % Compute adjacency
            adj = exp(-1*vecnorm(r(t_start_idx,:,bb)-r(t_start_idx,:,bbb),2,2).^2/(coeff_adj));

            % Update the adjacency matrix
            switch adj_type
                case 'spatial'
                    adj_trace{bb,bbb} = adj;
                    t_adj_trace{bb,bbb} = t_start_sec;
                    len_adj_trace{bb,bbb} = t_len_sec;
                case 'time'
                    adj_mat(bb,bbb) = sum(t_len_sec(adj>0.5)) / t_sum_sec;
            end
        end
    end
end
switch adj_type
    case 'spatial'
        varargout{1} = adj_trace;
        varargout{2} = t_adj_trace;
        varargout{3} = len_adj_trace;
    case 'time'
        varargout{1} = adj_mat;
end

