clc; clear; close all;

j = sqrt(-1);            % i or j

% Electromagnetic parameters
fc = 28e9;               % Carrier frequency (28GHz)
c = 3e8;                 % Speed of light (m/s)
lambda = c/fc;           % Wavelength (m/cycle)
BW = 50e3/2e8*fc;        % Bandwidth; only need a small bandwidth to see shift
fs = BW;                 % Sampling frequency (narrowband sampling < Nyquist)
k = 2*pi/lambda;         % wavenumber

% Simulation parameters 
sim_time = 1e-3;         % Simulation time (1 ms)
t = 0:1/fs:sim_time;     % Time vector
% split the signal into N parts. The blades will move by a little bit
% for each iteration.
N = 25;                  % num samples per update. 

% Radar parameters 
txPower = 10;            % Transmitted power (Watts)
noiseFigure = 6;         % Receiver noise figure (dB)

% Physical parameters 
L1 = 1;                  % length of the center of the rotor (m)
L2 = 2;                  % length of the blade + L1
rpm = 2e4;               % propeller rotations per min
rps = rpm/60*2*pi;       % rotation in rad per sec.
max_vel = rpm/60*2*pi*L2;% max velocity at rpm
R = 100;                 % Target range (100 m)
v = 0;                   % Target velocity (0 m/s)
N_scat = 50;             % N_scat scatters on each side of the origin for each
                         % propeller
N_prop = 4;              % number of propellers
N_blades = 2;            % number of blades
tx_pos = [0; 0; 0];      % Transmitter at origin
tx_vel = [0; 0; 0];      % Transmitter stationary


% contruct a uniformely distributed variable from L1 to L2 or only L2 if
% N_scat = 1
if N_scat > 1
    % bsp := blade scattering position
    bsp = L1:(L2-L1)/(N_scat-1):L2;
else
    bsp = L2;
end

% RCS init
rcs = zeros(1, N_blades*N_scat*N_prop);
rcs(1, :) = 0.05/N_scat;
% antenna is isotropic so angle doesn't matter. this is only used
% to indicate the number of reflectors.
scat_ang = zeros(1, N_blades*N_scat*N_prop);

% Full-system simulation using Phased Array
%-------------------------------------------------------------------------%
% Create Transmitter
tx = phased.Transmitter('PeakPower', txPower, 'Gain', 10);

% Create Antenna, Radiator and collector
ant = phased.IsotropicAntennaElement('BackBaffled',true);
txant = phased.Radiator('Sensor',ant,'PropagationSpeed',c,'OperatingFrequency',fc);
rxant = phased.Collector('Sensor',ant,'PropagationSpeed',c,'OperatingFrequency',fc);

% Create Target Model
target = phased.RadarTarget('MeanRCS', rcs, 'OperatingFrequency', fc, 'Model','Swerling2');

% Create Free Space Channel (Handles Delay and Doppler)
channel = phased.FreeSpace(...
    'SampleRate', fs, ...
    'TwoWayPropagation', true, ...
    'OperatingFrequency', fc);

% Create Receiver
rx = phased.ReceiverPreamp('Gain', 10, 'NoiseFigure', noiseFigure);

% Generate CW Signal
cw_waveform = exp(j*2*pi*fc*t).'; % Column vector signal
%-------------------------------------------------------------------------%

% Init data storage variables
num_iter = floor(length(t)/N); 
time_per_iter = sim_time/N;
rx_signal = zeros(size(cw_waveform));
theo_sig = zeros(size(cw_waveform));

% theoretical framework for pos and vel
% standardize notation: first underscore is subscript, second is superscript
x_u_o = [R;0;0]; 
x_t_o = [0;0;0];
v_t_o = [0;0;0];
R_u_o = eul2rotm([0 0 0]);
R_pi_u = eul2rotm([pi/2 pi/2 0]);
% the 0th and 2nd go counterclockwise and the other two, clockwise
prop_init_offset_angle = [pi/4;3*pi/4;5*pi/4;7*pi/4];
blade_init_offset_angle = [pi/4;3*pi/4;5*pi/4;7*pi/4];
y = 0.5; % vertical distance between the propellers and the center of the drone
x_bsp_o = zeros(3, N_prop*N_blades*N_scat);
v_bsp_o = zeros(3, N_prop*N_blades*N_scat);
mag_x_pi_u = 1;
vel_fac = 1;

% Main simulation loop
for i = 1:num_iter
    iter_range = ((i-1)*N+1:i*N);
	time = (i-1)*time_per_iter;
	idx = 1;
	
    % update scatterers' pos and vel
 %-------------------------------------------------------------------------% 
	% blade_ang is the angle of the first blade on the p-th propeller 
    blade_ang = blade_init_offset_angle + [time*rps;-time*rps;time*rps;-time*rps];
    
    if N_scat > 1
        for p = 1:N_prop
            theta = blade_ang(p);
            x_pi_u = [mag_x_pi_u*cos(prop_init_offset_angle(p));y;...
                mag_x_pi_u*sin(prop_init_offset_angle(p))];

            for b = 0:N_blades-1
                theta = theta + 2*pi*b/N_blades;

                for bsp = L1:(L2-L1)/(N_scat-1):L2
                    x_s_pi = bsp*[cos(theta);sin(theta);0];
                    x_s_u  = R_pi_u*x_s_pi+x_pi_u;
                    x_s_o  = R_u_o*x_s_u+x_u_o;

                    v_s_pi = vel_fac*bsp*[-sin(theta);cos(theta);0];
                    v_s_o = rps*R_u_o*R_pi_u*v_s_pi;

                    x_bsp_o(:,idx) = x_s_o;
                    v_bsp_o(:,idx) = v_s_o;
                    idx = idx + 1;
                end
            end
            vel_fac = -1*vel_fac;
        end
    else
        for p = 1:N_prop
            theta = blade_ang(p);
            x_pi_u = [mag_x_pi_u*cos(prop_init_offset_angle(p));y;...
                mag_x_pi_u*sin(prop_init_offset_angle(p))];

            for b = 0:N_blades-1
                theta = theta + 2*pi*b/N_blades;

                x_s_pi = bsp*[cos(theta);sin(theta);0];
                x_s_u  = R_pi_u*x_s_pi+x_pi_u;
                x_s_o  = R_u_o*x_s_u+x_u_o;

                v_s_pi = vel_fac*bsp*[-sin(theta);cos(theta);0];
                v_s_o = rps*R_u_o*R_pi_u*v_s_pi;

                x_bsp_o(:,idx) = x_s_o;
                v_bsp_o(:,idx) = v_s_o;
                idx = idx + 1;
            end
            vel_fac = -1*vel_fac;
        end
    end
 %-------------------------------------------------------------------------%

	% Phased Array package processing:
    % Transmit Signal
    tx_signal = cw_waveform(iter_range);
    
    % we need the number of scatters time this signal
    tx_signal = txant(tx_signal, scat_ang);
    prop_signal = channel(tx_signal, x_t_o, x_bsp_o, v_t_o, v_bsp_o);
    refl_signal = target(prop_signal, true);
    rx_signal(iter_range) = rxant(refl_signal, scat_ang);
end

% Receive the Signal
rx_signal = rx(rx_signal);

% theoretical:
 %-------------------------------------------------------------------------%
idx = 1;
vel_fac = 1;

for time=t
    blade_ang = blade_init_offset_angle + [time*rps;-time*rps;time*rps;-time*rps];
    res = 0;
    for p=1:N_prop
        x_pi_u = [mag_x_pi_u*cos(prop_init_offset_angle(p));y;...
                    mag_x_pi_u*sin(prop_init_offset_angle(p))];
        x_pi_o = R_u_o*x_pi_u + x_u_o;
        
        % sum over n 
        sn = 0;
        for b=1:N_blades
            omega=blade_ang(b);
            q = vel_fac*rps*R_u_o*R_pi_u*[-sin(omega);cos(omega);0];
            l = 4*pi*fc/c*(dot(q, x_pi_o/norm(x_pi_o)))*(time-2*norm(x_pi_o)/c);
            
            sn = sn + exp(1j*(L1+L2)*l/2)*sinc((L2-L1)*l/2);
        end

        res = res + sn*exp(-1j*4*pi*fc*norm(x_pi_o)/c);
        vel_fac = -1*vel_fac;
    end
    theo_sig(idx) = res;
    idx = idx + 1;
end
 %-------------------------------------------------------------------------%

% shift to baseband
rx_signal_bb = rx_signal .* exp(-1j*2*pi*fc*t).';
cw_waveform_bb = cw_waveform .* exp(-1j*2*pi*fc*t).';
theo_sig_bb = theo_sig .* exp(-1j*2*pi*fc*t).';

%plot time domain
%{
figure;
subplot(3,1,1);
plot(t, cw_waveform_bb, 'b');

subplot(3,1,2);
plot(t, rx_signal_bb, 'r');

subplot(3,1,3);
plot(t, theo_sig_bb, 'g');
%}


% Compute and plot spectrum
N = length(rx_signal_bb);

TX_FFT = abs(fftshift(fft(cw_waveform_bb, N)));
RX_FFT = abs(fftshift(fft(rx_signal_bb, N)));

TX_FFT = TX_FFT/max(TX_FFT);
RX_FFT = RX_FFT/max(RX_FFT);

% theoretical
rx_fft_theo = abs(fftshift(fft(theo_sig_bb, N)));
rx_fft_theo = rx_fft_theo/max(rx_fft_theo);

fft_spacing = -fs/2+1:fs/N:fs/2;

figure('Renderer', 'painters', 'Position', [10 10 900 600])
plot( fft_spacing/1e6, RX_FFT, 'r', fft_spacing/1e6, rx_fft_theo, 'g', fft_spacing/1e6, TX_FFT, 'b');
xline(-2*max_vel*fc/c/1e6);
xline(2*max_vel*fc/c/1e6);
legend('Rx simulation', 'Rx theoretical', 'Transmitted', 'Max + velocity', 'Max - velocity');
xlabel('Frequency (MHz)'); ylabel('Normalized Magnitude');
title('Baseband Frequency Spectrum');

