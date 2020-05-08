Fs = 8192;  % this is the sample rate
f_c = 1000;
msg_length = 13;
r = audiorecorder(Fs, 16, 1); % create an audiorecorder object
recordblocking(r, 3);     % record for 2 seconds

play(r)


y_r = double(getaudiodata(r, 'int16'));
t = [0:length(y_r)-1]*1/Fs;
plot(t, y_r);

