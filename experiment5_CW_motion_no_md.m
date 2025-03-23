clc; clear; close all;

%% Parameters
fc = 1e8;               % Carrier frequency (1 MHz)
c = 3e8;                 % Speed of light (m/s)
R = 200;                   % Target range (1000 m)
v = 0;                   % Target velocity (0 m/s)
lambda = c/fc;           % Wavelength
fs = 20*fc;              % Sampling frequency (at least Nyquist)
T = 4e-6;               % Simulation time
t = 0:1/fs:T;            % Time vector
txPower = 10;            % Transmitted power (Watts)
noiseFigure = 20;        % Receiver noise figure (dB)
k = 2*pi/lambda;         % wavenumber
j = sqrt(-1);            % i or j
L1 = 1;                  % length of the center of the rotor
L2 = 10;                  % length of the blade + L1
rps = 1/20;                % rotation in rad per sec.

N_scat = 10; % N_scat scatters on each side of the origin.
tx_pos = [0; 0; 0];        % Transmitter at origin
tx_vel = [0; 0; 0];        % Transmitter stationary
target_vel = zeros(3,2*N_scat);

R_vec = [R;0;0];
x1 = (N_scat:-1:1);
x2 = (1:N_scat);
m = (L2-L1)/(N_scat-1);
b = (L1*N_scat-L2)/(N_scat-1);  

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

% Propagate to Target and Reflect
% split the signal into N parts. The blades will move by a little bit
% for each iteration.
N = 20; % num samples per update
num_iter = length(t)/N; 
rx_signal = zeros(length(t), 1);
time_per_iter = T/N;

for i = 1:num_iter
    time = i*time_per_iter;
    blade_dir = [cos(time * rps);sin(time * rps);0];
    blade_dir = blade_dir/norm(blade_dir);
    test = pi/2-time * rps;
    target_pos = [R_vec-blade_dir*(m*x1+b) R_vec+blade_dir*(m*x2+b)];
    
    % Transmit Signal
    cw_sig_sub = zeros(length(t), 1);
    cw_sig_sub((i-1)*N+1:i*N) = cw_waveform((i-1)*N+1:i*N);
    tx_signal = tx(cw_sig_sub);
    
    % we need the number of scatters time this signal
    tx_signal = txant(tx_signal, scat_ang);

    %tx_sig_sub = zeros(1, ceil(1.5*length(t)));
    %tx_sig_sub((i-1)*N+1:i*N) = tx_signal((i-1)*N+1:i*N);
    prop_signal = channel(tx_signal, tx_pos, target_pos, tx_vel, target_vel);
    rx_signal = rx_signal + rxant(target(prop_signal), scat_ang);
end

% Receive the Signal
rx_signal = rx(rx_signal);

%% Plot Results
figure;

plot(t*1e6, real(rx_signal));
title('Received CW Signal - Exprimental:Blue, Theoretical:Red');
xlabel('Time (\mus)'); ylabel('Amplitude');
hold on

% theoretical
theo_res = zeros(1, length(t));

for i = 1:num_iter
    time = i*time_per_iter;
    gamma_0 = calc_gamma(k, pi/2-time * rps, R_vec);
    gamma_1 = calc_gamma(k, pi/2-time * rps + pi, R_vec);
    %{
    tsub = zeros(size(t));
    tsub((i-1)*N+1:i*N) = t((i-1)*N+1:i*N);
    %}
    r1 = zeros(size(theo_res));
    r2 = zeros(size(theo_res));
    r1((i-1)*N+1:i*N) = real(exp(j*2*pi*fc*(t((i-1)*N+1:i*N)-R/c)) * exp(j * (L1+L2) * gamma_0/2) * sinc((L2-L1)*gamma_0/2));
    r2((i-1)*N+1:i*N) = real(exp(j*2*pi*fc*(t((i-1)*N+1:i*N)-R/c)) * exp(j * (L1+L2) * gamma_1/2) * sinc((L2-L1)*gamma_1/2));
    theo_res = theo_res + (r1 + r2);
end

theo_res = 0.95 * theo_res * abs(max(rx_signal)) / abs (max(theo_res));

plot(t*1e6, theo_res, 'r');
xline(2*(R)/c * 1e6);

%{
r4 = 0.95 * real(exp(j*2*pi*fc*(t-R/c)));
r4 = r4 * abs(max(rx_signal)) / abs (max(r4));
plot(t*1e6, r4, 'g');
%}

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
    eul = [0 0 0];
    rot_po_o = eul2rotm(eul);
    
    p = sin(omega)*rot_po_o*[1;0;0] + cos(omega)*rot_po_o*[0;1;0];
    gamma = dot(-2*k*p, R/norm(R)); 
end



%}

