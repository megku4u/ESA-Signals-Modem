function start_idx = find_start_of_signal(x, sync)
% uses cross correlation based approach to find
% the start of the signal. it assumes that the
% signal in x has the sync signal at its very beginning


% The received signal includes a bunch of samples from before the
% transmission started so we need discard the samples from before
% the transmission started. 
% This block finds the start point of the signal based on
% the noise like signal that was transmitted 
% before the BPSK signal

% first make a coarse estimate of the start by looking for the first time
% index when the signal exceeds the rms value
signal_rms = rms(x);
coarse_idx = sort(find(x > signal_rms)); 

% back it up by some amount (we use the length of the sync here, 
% but it can be some other reasonable value) to ensure that the start is
% captured
coarse_idx = coarse_idx - length(sync);

[Ryx lags] = xcorr(x(1:coarse_idx+length(sync)*3), sync);
[mm ii] = max(abs(Ryx));
start_idx = lags(ii)+1;



end