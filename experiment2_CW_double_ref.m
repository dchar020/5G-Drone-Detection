clc; clear; close all;

%% Parameters
fc = 10e8;               % Carrier frequency (10 GHz)
c = 3e8;                 % Speed of light (m/s)
R = 15;                   % Target range (1000 m)
v = 0;                   % Target velocity (0 m/s)
lambda = c/fc;           % Wavelength
fs = 20*fc;              % Sampling frequency (at least Nyquist)
T = 8e-8;               % Simulation time
t = 0:1/fs:T;            % Time vector
rcs = [0.5 0.5];         % Radar cross-section (in m^2)
txPower = 10;            % Transmitted power (Watts)
noiseFigure = 20;        % Receiver noise figure (dB)
k = 2*pi/lambda;         % wavenumber
j = sqrt(-1);            % i or j
L1 = 1;                  % length of the center of the rotor
L2 = 2;                  % length of the blade + L1

tx_pos = [0; 0; 0];        % Transmitter at origin
target_pos = [R-L2 R+L2; 0 0; 0 0];    % Target at range R
tx_vel = [0; 0; 0];        % Transmitter stationary
target_vel = [0 0; 0 0; 0 0];    % Target moving along X

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
tx_signal = txant(tx_signal, [0 0]);

% Propagate to Target and Reflect
for i = 1:2
    prop_signal = channel(tx_signal, tx_pos, target_pos, tx_vel, target_vel);
    rx_signal = target(prop_signal);
    rx_signal = rxant(rx_signal, [0 0]);
end

% Receive the Signal
rx_signal = rx(rx_signal);

%% Plot Results
figure;

plot(t*1e6, real(rx_signal));
title('Received CW Signal - Exprimental:Blue, Theoretical:Red');
xlabel('Time (\mus)'); ylabel('Amplitude');
hold on

r1 = 1e1*real(exp(j*2*pi*fc*(t-R/c)) * sinc((L2)*k) * exp(j * (L1+L2) * k));
r2 = 1e1*real(exp(j*2*pi*fc*(t-R/c)) * sinc(-(L2)*k) * exp(-j * (L1+L2) * k));
r3 = (r1 + r2);

plot(t*1e6, r3*10, 'r');

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
