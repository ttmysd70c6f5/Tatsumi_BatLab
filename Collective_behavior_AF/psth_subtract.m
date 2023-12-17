function  ph = psth_subtract(times,times_shuffle,binsize, fs, ntrials, ntrials_shuffle, n_ite, triallen, smoothing, varargin)
% PSTH Computes the peri-stimulus time histogram from spike times.
% The routine plots the trial averaged spike rate as a function of time.
% H = PSTH(TIMES, BINSIZE, FS,NTRIALS,TRIALLEN)
% H = PSTH(TIMES, BINSIZE, FS,NTRIALS,TRIALLEN ,AXESHANDLE)
% TIMES - spike times (samples)
% BINSIZE - binwidth (ms)
% FS - sampling rate (hz)
% NTRIALS - number of trials
% TRIALLEN - length of a trial (samples)
% H - plot handle
%
% An example:
% %spike times can be specified in continuous time 
% %here we have 3 trials and a trial length of 1000 samples
% t = [10, 250, 900, 1300, 1600, 2405, 2900];
%
% %the same spike times can also be specified per trial
% t2 =[10, 250, 900, 300, 600, 405, 900];
% r = psth(t,10,1000,3,1000) ;
% r2 = psth(t2,10,1000,3,1000);
%
% Author: Rajiv Narayan
% askrajiv@gmail.com
% Boston University, Boston, MA
%
%
% <Additional feature>
% SMOOTHING: window size of movemean (samples)
% Added by Tatsumi Yoshida on 01/20/2023

h_color ='k';
nin=nargin;

error(nargchk(9,10, nin));

switch(nin)
 
 case 9 %no plot handle
  figure;
  h=gca;
  
 case 10
  if(ishandle(varargin{1}))
    h=varargin{1};
  else
    error('Invalid Plot handle');
  end

end

%Compute PSTH        
lastBin = binsize * ceil((triallen-1)*(1000/(fs*binsize))); % the last bin (msec)
edges = 0 : binsize : lastBin; % edge of psth (msec)
x_row = (mod(times-1,triallen)+1)*(1000/fs); % x axis of psth (msec). The first sample is set to 0-th bin
r = (histc(x_row,edges)*1000) / (ntrials*binsize); % Calculate firing rate (Fz). 1000 is multiplied to convert (1/msec) to (1/sec)

x_shuffle = (mod(times_shuffle-1,triallen)+1)*(1000/fs); % x axis of psth (msec). The first sample is set to 0-th bin
r_shuffle = (histc(x_shuffle,edges)*1000) / (n_ite*ntrials_shuffle*binsize); % Calculate firing rate (Fz). 1000 is multiplied to convert (1/msec) to (1/sec)

%Plot histogram        
axes(h);
ph=bar(edges(1:end-1),movmean(r(1:end-1)-r_shuffle(1:end-1),smoothing),'histc');
set(ph,'edgecolor',h_color,'facecolor',h_color);

