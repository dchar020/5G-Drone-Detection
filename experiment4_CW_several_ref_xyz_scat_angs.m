clc; clear; close all;

%% Parameters
fc = 1e8;               % Carrier frequency (1 MHz)
c = 3e8;                 % Speed of light (m/s)
R = 1000;                   % Target range (1000 m)
v = 0;                   % Target velocity (0 m/s)
lambda = c/fc;           % Wavelength
fs = 20*fc;              % Sampling frequency (at least Nyquist)
T = 15e-6;               % Simulation time
t = 0:1/fs:T;            % Time vector
txPower = 10;            % Transmitted power (Watts)
noiseFigure = 20;        % Receiver noise figure (dB)
k = 2*pi/lambda;         % wavenumber
j = sqrt(-1);            % i or j
L1 = 1;                  % length of the center of the rotor
L2 = 10;                  % length of the blade + L1

N_scat = 10; % N_scat scatters on each side of the origin.
tx_pos = [0; 0; 0];        % Transmitter at origin
tx_vel = [0; 0; 0];        % Transmitter stationary

blade_dir = [1;1;1];
blade_dir = blade_dir/norm(blade_dir);
R_vec = [R;0;0];
x1 = (N_scat:-1:1);
x2 = (1:N_scat);

m = (L2-L1)/(N_scat-1);
b = (L1*N_scat-L2)/(N_scat-1);

target_pos = [R_vec-blade_dir*(m*x1+b) R_vec+blade_dir*(m*x2+b)];


%target_pos = [R_vec-blade_dir R_vec+blade_dir];
target_vel = zeros(3,2*N_scat);    % Target moving along X

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

% we need the number of scatters time this signal
tx_signal = txant(tx_signal, scat_ang);

% Propagate to Target and Reflect
prop_signal = channel(tx_signal, tx_pos, target_pos, tx_vel, target_vel);
rx_signal = target(prop_signal, true);
rx_signal = rxant(rx_signal, scat_ang);

% Receive the Signal
rx_signal = rx(rx_signal);

%% Plot Results
figure;

start_window = 15000;
end_window = 15200;

plot(t(start_window:end_window)*1e6, real(rx_signal(start_window:end_window)));
title('Received CW Signal - Exprimental:Blue, Theoretical:Red');
xlabel('Time (\mus)'); ylabel('Amplitude');
hold on

% theoretical
gamma_0 = calc_gamma(k, pi/4, R_vec);
gamma_1 = calc_gamma(k, 5*pi/4, R_vec);

r1 = real(exp(j*2*pi*fc*(t-R/c)) * exp(j * (L1+L2) * gamma_0/2) * sinc((L2-L1)/2*gamma_0));
r2 = real(exp(j*2*pi*fc*(t-R/c)) * exp(j * (L1+L2) * gamma_1/2) * sinc((L2-L1)/2*gamma_1));
r3 = (r1 + r2);

r3 = 0.95 * r3 * abs(max(rx_signal)) / abs (max(r3));

plot(t(start_window:end_window)*1e6, r3(start_window:end_window), 'r');
%xline(2*(R)/c * 1e6);

r4 = 0.95 * real(exp(j*2*pi*fc*(t-R/c)));
r4 = r4 * abs(max(rx_signal)) / abs (max(r4));
plot(t(start_window:end_window)*1e6, r4(start_window:end_window), 'g');

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


function gamma = calc_gamma(k, omega, R)
    eul = [0 -pi/4 0];
    rot_po_o = eul2rotm(eul);
    
    p = sin(omega)*rot_po_o*[1;0;0] + cos(omega)*rot_po_o*[0;1;0];
    gamma = dot(-2*k*p, R/norm(R)); 
end



%}

