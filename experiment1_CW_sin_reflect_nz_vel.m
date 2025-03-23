clc; clear; close all;

%% Parameters
fc = 10e9;               % Carrier frequency (10 GHz)
c = 3e8;                 % Speed of light (m/s)
R = 100;                % Target range (1000 m)
v = 50;                  % Target velocity (50 m/s)
lambda = c/fc;           % Wavelength
BW = 10e3;               % only need 10 kHz bandwidth
fs = BW;               % Sampling frequency (at least Nyquist)
T = 1e-3;                % Simulation time
t = 0:1/fs:T;            % Time vector
rcs = 1;                 % Radar cross-section (in m^2)
txPower = 1;             % Transmitted power (Watts)
noiseFigure = 5;         % Receiver noise figure (dB)

tx_pos = [0; 0; 0];        % Transmitter at origin
target_pos = [R; 0; 0];    % Target at range R
tx_vel = [0; 0; 0];        % Transmitter stationary
target_vel = [v; 0; 0];    % Target moving along X

%% Create Transmitter
tx = phased.Transmitter('PeakPower', txPower, 'Gain', 40);

%% Create Target Model
target = phased.RadarTarget('MeanRCS', rcs, 'OperatingFrequency', fc);

%% Create Free Space Channel (Handles Delay and Doppler)
channel = phased.FreeSpace(...
    'SampleRate', fs, ...
    'TwoWayPropagation', true, ...
    'OperatingFrequency', fc);

%% Create Receiver
rx = phased.ReceiverPreamp('Gain', 40, 'NoiseFigure', noiseFigure);

%% Generate CW Signal
cw_waveform = exp(1j*2*pi*fc*t).';  % Column vector signal

% Transmit Signal
tx_signal = tx(cw_waveform);

% Propagate to Target and Reflect
prop_signal = channel(tx_signal, tx_pos, target_pos, tx_vel, target_vel);
rx_signal = target(prop_signal);

% Receive the Signal
rx_signal = rx(rx_signal);

% shift the signal to baseband; element wise operation
rx_signal_bb = rx_signal .* exp(-1j*2*pi*fc*t).';
tx_signal_bb = tx_signal .* exp(-1j*2*pi*fc*t).';

hold on
rx_signal_bb = rx_signal_bb/abs(max(rx_signal_bb));
tx_signal_bb = tx_signal_bb/abs(max(tx_signal_bb));

plot(t, rx_signal_bb);
plot(t, tx_signal_bb, 'r');
hold off

% take the fft
N = 128;
TX_FFT = abs(fftshift(fft(tx_signal_bb, N)));
RX_FFT = abs(fftshift(fft(rx_signal_bb, N)));

TX_FFT = TX_FFT/max(TX_FFT);
RX_FFT = RX_FFT/max(RX_FFT);

figure;
plot(-fs/2+1:fs/N:fs/2, TX_FFT, 'b', -fs/2+1:fs/N:fs/2, RX_FFT, 'r');
legend('Transmitted', 'Received');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Baseband Frequency Spectrum (Doppler Shift Observed)');
xline(-2*v*fc/c);

%% Plot Results
%{
figure;
subplot(2,1,1);
plot(t*1e6, real(tx_signal));
title('Transmitted CW Signal');
xlabel('Time (\mus)'); ylabel('Amplitude');

subplot(2,1,2);
plot(t*1e6, real(rx_signal));
title('Received CW Signal (with Delay and Doppler Shift)');
xlabel('Time (\mus)'); ylabel('Amplitude');
%}

% Compute and plot spectrum
%{
N = length(t);
f = linspace(-fs/2, fs/2, N);
TX_FFT = abs(fftshift(fft(tx_signal, N)));
RX_FFT = abs(fftshift(fft(rx_signal, N)));

TX_FFT = TX_FFT/max(TX_FFT);
RX_FFT = RX_FFT/max(RX_FFT);

%[max, argmax_tx] = max(TX_FFT(floor(N/2)+1:N));
[max, argmax_rx] = max(RX_FFT(floor(N/2)+1:N));

figure;
plot(f/1e6, TX_FFT, 'b', f/1e6, RX_FFT, 'r');
legend('Transmitted', 'Received');
xlabel('Frequency (MHz)'); ylabel('Magnitude');
title('Frequency Spectrum (Doppler Shift Observed)');

fd = 2*v*fc/c*1e-6;
xline(10000-fd);
xlim([9920, 10080]);
ylim([0 0.01]);
%}