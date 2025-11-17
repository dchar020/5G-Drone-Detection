clc; clear; close all;

%% Parameters
fc = 1e8;               % Carrier frequency (1 GHz)
c = 3e8;                 % Speed of light (m/s)
R = 15;                   % Target range (1000 m)
v = 0;                   % Target velocity (0 m/s)
lambda = c/fc;           % Wavelength
fs = 20*fc;              % Sampling frequency (at least Nyquist)
T = 30e-8;               % Simulation time
t = 0:1/fs:T;            % Time vector
txPower = 10;            % Transmitted power (Watts)
noiseFigure = 20;        % Receiver noise figure (dB)
k = 2*pi/lambda;         % wavenumber
j = sqrt(-1);            % i or j
L1 = 1;                  % length of the center of the rotor
L2 = 2;                  % length of the blade + L1

N_scat = 30; % N_scat scatters on each side of the origin.
tx_pos = [0; 0; 0];        % Transmitter at origin
target_pos = [R-(L2-L1)/N_scat*(N_scat:-1:1) R+(L2-L1)/N_scat*(1:N_scat); ...
    zeros(1,2*N_scat); zeros(1, 2*N_scat)];
tx_vel = [0; 0; 0];        % Transmitter stationary
target_vel = [zeros(3,2*N_scat)];    % Target moving along X

rcs = zeros(1, N_scat * 2);
rcs(1, :) = 0.5;
scat_ang = zeros(1, N_scat * 2);

%% Create Transmitter
tx = phased.Transmitter('PeakPower', txPower, 'Gain', 40);

%% Create Antenna, Radiator and collector
ant = phased.IsotropicAntennaElement('BackBaffled',true);
txant = phased.Radiator('Sensor',ant,'PropagationSpeed',c,'OperatingFrequency',fc);
rxant = phased.Collector('Sensor',ant,'PropagationSpeed',c,'OperatingFrequency',fc);

%% Create Target Model
target = phased.RadarTarget('MeanRCS', rcs, 'OperatingFrequency', fc, 'Model','Nonfluctuating');

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

% we need the number of scatters time this signal
tx_signal = txant(tx_signal, scat_ang);

% Propagate to Target and Reflect
prop_signal = channel(tx_signal, tx_pos, target_pos, tx_vel, target_vel);
rx_signal = target(prop_signal);
rx_signal = rxant(rx_signal, scat_ang);

% Receive the Signal
rx_signal = rx(rx_signal);

%% Plot Results
figure;

plot(t*1e6, real(rx_signal));
title('Received CW Signal - Exprimental:Blue, Theoretical:Red');
xlabel('Time (\mus)'); ylabel('Amplitude');
hold on

r1 = 1e1*real(exp(j*2*pi*fc*(t-R/c)) * sinc((L2+L1)*k) * exp(j * (L1+L2) * k));
r2 = 1e1*real(exp(j*2*pi*fc*(t-R/c)) * sinc(-(L2+L1)*k) * exp(-j * (L1+L2) * k));
r3 = (r1 + r2);

r3 = - r3 / max(abs(r3)) * max(rx_signal);

plot(t*1e6, r3, 'r');
xline(2*(R+L2)/c * 1e6);

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
%}
