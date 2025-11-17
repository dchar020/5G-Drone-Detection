clc; clear; close all;

%% Parameters
fc = 10e8;               % Carrier frequency (10 GHz)
c = 3e8;                 % Speed of light (m/s)
R = 1;                   % Target range (1000 m)
v = 0;                   % Target velocity (0 m/s)
lambda = c/fc;           % Wavelength
bandwidth = 10e6;
fs = bandwidth;              % Sampling frequency (at least Nyquist)
T = 1e-1;               % Simulation time
t = 0:1/fs:T;            % Time vector
rcs = 1;                 % Radar cross-section (in m^2)
txPower = 10;            % Transmitted power (Watts)
noiseFigure = 20;        % Receiver noise figure (dB)
k = 2*pi/lambda;
j = sqrt(-1);

tx_pos = [0; 0; 0];        % Transmitter at origin
target_pos = [R; 0; 0];    % Target at range R
tx_vel = [0; 0; 0];        % Transmitter stationary
target_vel = [v; 0; 0];    % Target moving along X

%% Create Transmitter
tx = phased.Transmitter('PeakPower', txPower, 'Gain', 40);

%% Create Target Model
target = phased.RadarTarget('MeanRCS', rcs, 'OperatingFrequency', fc, 'Model','Swerling1');

%% Create Free Space Channel (Handles Delay and Doppler)
channel = phased.FreeSpace(...
    'SampleRate', fs, ...
    'TwoWayPropagation', true, ...
    'OperatingFrequency', fc);

%% Create Receiver
rx = phased.ReceiverPreamp('Gain', 40, 'NoiseFigure', noiseFigure);

%% Generate CW Signal
cw_waveform = exp(j*2*pi*fc*t).'; % Column vector signal

% Transmit Signal
tx_signal = tx(cw_waveform);

% Propagate to Target and Reflect
prop_signal = channel(tx_signal, tx_pos, target_pos, tx_vel, target_vel);
rx_signal = target(prop_signal, true);

% Receive the Signal
rx_signal = rx(rx_signal);

%% Plot Results
figure;
subplot(2,1,1);
plot(t*1e6, real(tx_signal));
title('Transmitted CW Signal');
xlabel('Time (\mus)'); ylabel('Amplitude');

subplot(2,1,2);
plot(t*1e6, real(rx_signal));
title('Received CW Signal - Exprimental:Blue, Theoretical:Red');
xlabel('Time (\mus)'); ylabel('Amplitude');
hold on

plot(t*1e6, 100*real(exp(j*2*pi*fc*(t-R/c))), 'r');
xline(2*R/c * 1e6,'g')

hold off


% Compute and plot spectrum
%{
N = length(t);
f = linspace(-fs/2, fs/2, N);
TX_FFT = abs(fftshift(fft(tx_signal, N)));
RX_FFT = abs(fftshift(fft(rx_signal, N)));

figure;
plot(f/1e6, TX_FFT, 'b', f/1e6, RX_FFT, 'r');
legend('Transmitted', 'Received');
xlabel('Frequency (MHz)'); ylabel('Magnitude');
title('Frequency Spectrum (Doppler Shift Observed)');
%}
