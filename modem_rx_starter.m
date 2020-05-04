
%load short_modem_rx.mat
load long_modem_rx.mat

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

plot(y_t);
% create time vector
t = linspace(0, size(y_t,1)-1, size(y_t,1))';
t = t./Fs;

figure()
DTFT_fftbased(y_t);

% Scale down by factor of G
y_t = y_t./500;

%multiply cosine
y_cos = y_t.*cos(2.*pi.*f_c.*t);

figure()
DTFT_fftbased(y_cos);

% Apply low pass filter
W = 480;
h = W/pi*sinc(W/pi*t);

y_read = conv(y_cos, h);
figure()
DTFT_fftbased(y_read);

figure()
plot(y_read);

%% 
% Approximate signal
y = y_read(50:100:end);
y = y(1:msg_length*8);
y(y > 0) = 1;
y(y < 0) = 0;

figure()
plot(y);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_d = y;
% convert to a string assuming that x_d is a vector of 1s and 0s
% representing the decoded bits
BitsToString(x_d)

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
