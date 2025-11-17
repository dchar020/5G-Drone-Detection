clc; clear; close all;

%% Parameters
fc = 1e9;               % Carrier frequency (1 GHz)
c = 3e8;                 % Speed of light (m/s)
R = 100;                   % Target range (1000 m)
v = 0;                   % Target velocity (0 m/s)
lambda = c/fc;           % Wavelength
fs = 2*fc;              % Sampling frequency (must be at least Nyquist
                      % for the phased package to handle the signal
                      % properly.
sim_time = 1e-6;               % Simulation time
% i don't care what happens before the wave reflects so 
t = 0:1/fs:sim_time;            % Time vector
txPower = 10;            % Transmitted power (Watts)
noiseFigure = 20;        % Receiver noise figure (dB)
k = 2*pi/lambda;         % wavenumber
j = sqrt(-1);            % i or j
L1 = 1;                  % length of the center of the rotor
L2 = 10;                  % length of the blade + L1
rps = 20;                % rotation in rad per sec.

N_scat = 2; % N_scat scatters on each side of the origin.
tx_pos = [0; 0; 0];        % Transmitter at origin
tx_vel = [0; 0; 0];        % Transmitter stationary
target_vel = zeros(3,2*N_scat);

R_vec = [R;0;0];
x1 = (N_scat:-1:1);
x2 = (1:N_scat);
m = (L2-L1)/(N_scat-1);
b = (L1*N_scat-L2)/(N_scat-1);  
y1 = (-1:1/N_scat:-L1/L2);
y2 = (L1/L2:1/N_scat:1);

rcs = zeros(1, N_scat * 2);
rcs(1, :) = 0.001;
% antenna is isotropic so angle doesn't matter. this is only used
% to indicate the number of reflectors.
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

blade_dir = [1;0;0];
target_pos = [R_vec-blade_dir*(m*x1+b) R_vec+blade_dir*(m*x2+b)];
rx_signal = zeros(size(cw_waveform)); 

tic
for i = 1:length(cw_waveform)
    tx_signal = tx(cw_waveform(i)); 
    tx_signal = txant(tx_signal, scat_ang);
    prop_signal = channel(tx_signal, tx_pos, target_pos, tx_vel, target_vel);
    refl_sig = target(prop_signal);
    rx_signal(i) = rxant(refl_sig, scat_ang);
end
toc

% Receive the Signal
rx_signal = rx(rx_signal);

%% Plot Results
figure;

plot(t*1e6, real(rx_signal));
title('Received CW Signal - Exprimental:Blue, Theoretical:Red');
xlabel('Time (\mus)'); ylabel('Amplitude');
hold on
xline(2*(R)/c * 1e6);
hold off


