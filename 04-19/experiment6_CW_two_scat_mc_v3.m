clc; clear; close all;

%% Parameters
rng(20250316);
fc = 1e9;               % Carrier frequency (1 GHz)
c = 3e8;                 % Speed of light (m/s)
R = 100;                   % Target range (100 m)
v = 0;                   % Target velocity (0 m/s)
lambda = c/fc;           % Wavelength
BW = 100e3;               % only need a small bandwidth to see shift
fs = BW;              % Sampling frequency
                      % for the phased package to handle the signal
                      % properly.
%sim_time = 1e-6;               % Simulation time
sim_time = 1e-3;
% i don't care what happens before the wave reflects so 
t = 0:1/fs:sim_time;            % Time vector
txPower = 10;            % Transmitted power (Watts)
noiseFigure = 0;        % Receiver noise figure (dB)
k = 2*pi/lambda;         % wavenumber
j = sqrt(-1);            % i or j
L1 = 1;                  % length of the center of the rotor
L2 = 2;                  % length of the blade + L1
rps = 2e4/60*2*pi;                % rotation in rad per sec.
max_vel = rps*L2;

multiple_scat = true;
N_scat = 100; % N_scat scatters on each side of the origin.
tx_pos = [0; 0; 0];        % Transmitter at origin
tx_vel = [0; 0; 0];        % Transmitter stationary
target_vel = zeros(3,2*N_scat);

R_vec = [R;0;0];

if multiple_scat
    x1 = (N_scat:-1:1);
    x2 = (1:N_scat);
    m = (L2-L1)/(N_scat-1);
    b = (L1*N_scat-L2)/(N_scat-1);  
    y1 = (-L2:(L2-L1)/(N_scat-1):-L1)/L2;
    y2 = (L1:(L2-L1)/(N_scat-1):L2)/L2;
else
    x1 = 0;
    x2 = 0;
    m = 0; 
    b = L2;
    y1 = -1;
    y2 = 1;
end

rcs = zeros(1, N_scat * 2);
rcs(1, :) = 0.1;
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

% Propagate to Target and Reflect
% split the signal into N parts. The blades will move by a little bit
% for each iteration.
N = 25; % num samples per update. LEAVE THIS AT 25
num_iter = floor(length(t)/N); 
time_per_iter = sim_time/N;
rx_signal = zeros(size(cw_waveform));
theo_sig = zeros(size(cw_waveform));

% debug
tar_pos_x_store = size((t));
tar_pos_y_store = size((t));
tar_vel_x_store = size((t));
tar_vel_y_store = size((t));

for i = 1:num_iter
    iter_range = ((i-1)*N+1:i*N);
    % update scatterers
    time = (i)*time_per_iter;
    blade_dir = [cos(time*rps);sin(time*rps);0];
    blade_dir = blade_dir/norm(blade_dir);
    target_pos = [R_vec-blade_dir*(m*x1+b) R_vec+blade_dir*(m*x2+b)];
    tg_dir = [-sin(time*rps); cos(time*rps); 0];
    target_vel = [L2*rps*tg_dir*y1 L2*rps*tg_dir*y2];

    %{
    tar_pos_x_store(iter_range) = target_pos(1,1);
    tar_pos_y_store(iter_range) = target_pos(2,1);
    tar_vel_x_store(iter_range) = target_vel(1,1);
    tar_vel_y_store(iter_range) = target_vel(2,1);
    %}

    % Transmit Signal
    tx_signal = cw_waveform((i-1)*N+1:i*N);
    
    % we need the number of scatters time this signal
    tx_signal = txant(tx_signal, scat_ang);
    prop_signal = channel(tx_signal, tx_pos, target_pos, tx_vel, target_vel);
    refl_signal = target(prop_signal);
    rx_signal((i-1)*N+1:i*N) = rxant(refl_signal, scat_ang);
end

%{
hold on
plot(1:length(t), tar_pos_x_store);
plot(1:length(t), tar_pos_y_store);
hold off
%}

% theoretical
idx = 1;
for time=t
    omega1 = time*rps;
    omega2 = time*rps+pi;
    q1 = rps*[-sin(omega1); cos(omega1); 0];
    q2 = rps*[-sin(omega2); cos(omega2); 0];
    l1 = 4*pi*fc/c*(dot(q1, R_vec/R))*(time-2*R/c);
    l2 = 4*pi*fc/c*(dot(q2, R_vec/R))*(time-2*R/c);
    theo_sig(idx) = exp(j*2*pi*fc*(time-2*R/c)) * ( ...
        exp(j*(L1+L2)*l1/2)*sinc((L2-L1)*l1/2) + ...
        exp(j*(L1+L2)*l2/2)*sinc((L2-L1)*l2/2) );
    idx = idx + 1;
end

% Receive the Signal
rx_signal = rx(rx_signal);

% Compute and plot spectrum
start_idx = ceil(2*(R+L2)/c * 1e6 * length(t));
end_idx = length(t);
N_samples = end_idx - start_idx;

% shift to baseband
rx_signal_bb = rx_signal .* exp(-1j*2*pi*fc*t).';
cw_waveform_bb = cw_waveform .* exp(-1j*2*pi*fc*t).';
theo_sig_bb = theo_sig .* exp(-1j*2*pi*fc*t).';

%plot(t,rx_signal_bb);

%cw_waveform = cw_waveform / abs(max(cw_waveform));
%rx_signal = rx_signal / abs(max(rx_signal));

N = 2^10;

TX_FFT = abs(fftshift(fft(cw_waveform_bb, N)));
RX_FFT = abs(fftshift(fft(rx_signal_bb, N)));

TX_FFT = TX_FFT/max(TX_FFT);
RX_FFT = RX_FFT/max(RX_FFT);

% theoretical
rx_fft_theo = abs(fftshift(fft(theo_sig_bb, N)));
rx_fft_theo = rx_fft_theo/max(rx_fft_theo);

fft_spacing = -fs/2+1:fs/N:fs/2;

figure;
%plot(f/1e6, TX_FFT, 'b');
plot(fft_spacing, TX_FFT, 'b', fft_spacing, RX_FFT, 'r', fft_spacing, rx_fft_theo, 'g');
xline(-2*max_vel*fc/c);
xline(2*max_vel*fc/c);
legend('Transmitted', 'Received', 'Theoretical', 'Max + velocity', 'Max - velocity');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Frequency Spectrum (Doppler Shift Observed)');


hold on

% show STFT with window = 7 us
%{
figure;
stft(rx_signal,fs, FFTLength=2^12);
%spectrogram(rx_signal);
%}

hold off
%xlabel('Frequency (MHz)'); ylabel('Magnitude');
%title('Theoretical Frequency Spectrum (Doppler Shift Observed)');
%xlim([900 1100]);
