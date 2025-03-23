clc; clear; close all;

%% Parameters
fc = 1e9;               % Carrier frequency (1 GHz)
c = 3e8;                 % Speed of light (m/s)
R = 100;                   % Target range (1000 m)
v = 0;                   % Target velocity (0 m/s)
lambda = c/fc;           % Wavelength
fs = 6*fc;              % Sampling frequency (must be at least Nyquist
                      % for the phased package to handle the signal
                      % properly.
%sim_time = 1e-6;               % Simulation time
sim_time = 1e-5;
% i don't care what happens before the wave reflects so 
t = 0:1/fs:sim_time;            % Time vector
txPower = 10;            % Transmitted power (Watts)
noiseFigure = 20;        % Receiver noise figure (dB)
k = 2*pi/lambda;         % wavenumber
j = sqrt(-1);            % i or j
L1 = 1;                  % length of the center of the rotor
L2 = 10;                  % length of the blade + L1
rps = 2e4/60*2*pi;                % rotation in rad per sec.

multiple_scat = false;
N_scat = 1; % N_scat scatters on each side of the origin.
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
    y1 = -L2;
    y2 = L2;
end

rcs = zeros(1, N_scat * 2);
rcs(1, :) = 2;
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
target = phased.RadarTarget('MeanRCS', rcs, 'OperatingFrequency', fc, 'Model','Swerling3');

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
N = 20; % num samples per update
num_iter = length(t)/N; 
time_per_iter = sim_time/N;
%tx_signal = cw_waveform(1:200);
rx_signal = zeros(size(cw_waveform));
theo_sig = zeros(size(cw_waveform));

for i = 1:num_iter
    % update scatterers
    time = i*time_per_iter;
    blade_dir = [cos(time*rps);sin(time*rps);0];
    blade_dir = blade_dir/norm(blade_dir);
    target_pos = [R_vec-blade_dir*(m*x1+b) R_vec+blade_dir*(m*x2+b)];
    tg_dir = [sin(time*rps); cos(time*rps); 0];
    target_vel = [(L2-L1)*rps*tg_dir*y1 (L2-L1)*rps*tg_dir*y2];
    
    % Transmit Signal
    tx_signal = cw_waveform((i-1)*N+1:i*N);
    
    % we need the number of scatters time this signal
    tx_signal = txant(tx_signal, scat_ang);
    prop_signal = channel(tx_signal, tx_pos, target_pos, tx_vel, target_vel);
    refl_signal = target(prop_signal, true);
    rx_signal((i-1)*N+1:i*N) = rxant(refl_signal, scat_ang);

    % theoretical
    time_theo = t((i-1)*N+1:i*N);
    omega = [time*rps;time*rps+pi];
    q1 = rps*[-sin(omega(1)); cos(omega(1)); 0];
    q2 = rps*[-sin(omega(2)); cos(omega(2)); 0];
    l = [4*pi*fc/c*(dot(q1, R_vec/R))*(time-2*R/c), ...
         4*pi*fc/c*(dot(q2, R_vec/R))*(time-2*R/c)];
    theo_sig((i-1)*N+1:i*N) = exp(j*2*pi*fc*(time_theo-2*R/c)) * ( ...
        exp(j*(L1+L2)*l(1)/2)*sinc((L2-L1)*l(1)/2) + ...
        exp(j*(L1+L2)*l(2)/2)*sinc((L2-L1)*l(2)/2) );
end

% Receive the Signal
rx_signal = rx(rx_signal);

%% Plot Results
figure;

plot(t*1e6, real(rx_signal));
title('Received CW Signal - Exprimental:Blue, Theoretical:Red');
xlabel('Time (\mus)'); ylabel('Amplitude');
hold on
xline(2*(R+L2)/c * 1e6);
hold off 

% Compute and plot spectrum
start_idx = ceil(2*(R+L2)/c * 1e6 * length(t));
end_idx = length(t);
N = end_idx - start_idx;
f = linspace(-fs/2, fs/2, N);
cw_waveform = cw_waveform / abs(max(cw_waveform));
rx_signal = rx_signal / abs(max(rx_signal));

TX_FFT = abs(fftshift(fft(cw_waveform(start_idx:end_idx), N)));
RX_FFT = abs(fftshift(fft(rx_signal(start_idx:end_idx), N)));

TX_FFT = TX_FFT/max(TX_FFT);
RX_FFT = RX_FFT/max(RX_FFT);

figure;
%plot(f/1e6, TX_FFT, 'b');
plot(f/1e6, TX_FFT, 'b', f/1e6, RX_FFT, 'r');
legend('Transmitted', 'Received');
xlabel('Frequency (MHz)'); ylabel('Magnitude');
title('Frequency Spectrum (Doppler Shift Observed)');

xlim([900 1100]);

% show STFT with window = 7 us
figure;
stft(rx_signal,fs, FFTLength=2^12);
%spectrogram(rx_signal);

% theoretical
rx_fft_theo = abs(fftshift(fft(theo_sig(start_idx:end_idx), N)));
rx_fft_theo = rx_fft_theo/max(rx_fft_theo);

figure;
plot(f/1e6, rx_fft_theo);
xlabel('Frequency (MHz)'); ylabel('Magnitude');
title('Theoretical Frequency Spectrum (Doppler Shift Observed)');
xlim([900 1100]);
