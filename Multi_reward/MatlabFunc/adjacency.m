function adj = adjacency(x,th_adj)
% x: 3d distance [t,d]
% th_adj: threshold for adjacency. adjacency equals to 0.5 if the distance
% equals to th_adj
%
% OUTPUT
% adj: adjacency rangin between [0,1]. adj = 1 means two data is at the
% same position.

% th_adj = 0.5;
coeff_adj = th_adj^2/(log(2)); % scaling coefficient
adj = exp(-1*vecnorm(x,2,2).^2/(coeff_adj)); % x: [t,d], adjacency: [t,]

