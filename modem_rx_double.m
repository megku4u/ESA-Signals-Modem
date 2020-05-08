close all

% Choose which message to decode
%load short_modem_rx.mat
load double_maybe_rx.mat

% The received signal includes a bunch of samples from before the
% transmission started so we need discard the samples from before
% the transmission started. 

start_idx = find_start_of_signal(y_r,x_sync);

% start_idx now contains the location in y_r where x_sync begins
% we need to offset by the length of x_sync to only include the signal
% we are interested in

y_t = y_r(start_idx+length(x_sync):end); % y_t is the signal which starts at the beginning of the transmission


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% create time vector
t = linspace(0, size(y_t,1)-1, size(y_t,1))';
t = t./Fs;

t_sinc = linspace(0,2,size(y_t,1));

% figure
% plot(y_t);
% title("Raw Received Signal")
% 
% figure
% title("DTFT of Raw Signal")
% DTFT_fftbased(y_t);

%%
% Scale down by factor of G
y_t = y_t./(4*10^6);

%multiply by cosine
y_cos = y_t.*cos(2.*pi.*f_c.*t);
y_sin = y_t.*sin(2.*pi.*f_c.*t);


% figure
% title("DTFT of Signal Multiplied by Cosine")
% DTFT_fftbased(y_cos);
% 
% figure
% title("DTFT of Signal Multiplied by Sine")
% DTFT_fftbased(y_sin);


%%
% Apply low pass filter
W = f_c/2;
h = W/pi*sinc(W/pi*t_sinc);

% Convolve the signal with the low-pass filter and plot the output
y_cos_read = conv(y_cos, h);
y_sin_read = -conv(y_sin, h);


figure
title("DTFT of Cos Convolved Signal")
DTFT_fftbased(y_cos_read);

figure
title("DTFT of Sin Convolved Signal")
DTFT_fftbased(y_sin_read);

figure(6)
plot(y_cos_read(1:100*msg_length*8));
title("Filtered Signal")

%% 
% Approximate signal - take out the noise and the extra length at the end.
% y is the vector of 1s and 0s that represent the message in bits.
y_cos = y_cos_read(50:100:end);
y_sin = y_sin_read(50:100:end);

y_cos = y_cos(1:msg_length*8);
y_sin = y_sin(1:msg_length*8);

y_cos(y_cos > 0) = 1;
y_cos(y_cos < 0) = 0;

y_sin(y_sin > 0) = 1;
y_sin(y_sin < 0) = 0;

% figure(7)
% plot(y);
% title("Decoded Signal")

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert to a string assuming that x_d is a vector of 1s and 0s
% representing the decoded bits
BitsToString(y_cos)
BitsToString(y_sin)

%% FUNCTIONS
function [X, O] = DTFT_fftbased(x)
    % plots the magnitude of the Fourier transform of the signal x
    % which is assumed to originate from a Continous-time signal 
    % sampled with frequency fs
    % the function returns X and f.
    % In other words, this function plots the FT of the DT signal x
    % with the frequency axis labeled as if it were the original CT signal
    % 
    % X contains the frequency response
    % O contains the frequency samples 

    N = length(x);

    X = fftshift(fft(x));
    O = linspace(-pi, pi - 2*pi/length(x), length(x));
    subplot(211);
    plot(O, abs(X));
    xlabel('\Omega');
    ylabel('|X(j\Omega)|');
    subplot(212);
    plot(O, unwrap(angle(X)));
    xlabel('\Omega');
    ylabel('\angle X(j\Omega)');
    
end
