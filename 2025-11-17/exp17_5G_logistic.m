clc; clear; close all;

%% Parameters
%rng(20250323);
fc = 50e9;                % Carrier frequency (20GHz)
c = 3e8;                 % Speed of light (m/s)
R = 100;                 % Target range (100 m)
v = 0;                   % Target velocity (0 m/s)
lambda = c/fc;           % Wavelength

% OFDM Parameters 
T_symbol = 1/(960e3); % OFDM symbol duration for 5G. Smallest possible value
N_s = 14;         % num subcarriers
SCS = 1/T_symbol; % subcarrier spacing
BW = SCS*(N_s+1);
fs = BW;                 % Sampling frequency (narrowband sampling < Nyquist)

pulse_period = 2e-6; %2e-6
PRF = 1/pulse_period;
N = round(fs/PRF); % num samples per update. 
num_iter = 2e3; % number of pulses

t_symbol = 0:1/fs:T_symbol; % interval for the first symbol
 
txPower = 3.16;            % Transmitted power (Watts)
noiseFigure = 20;        % Receiver noise figure (dB)
G_rx = 40;               % Receiver gain (dB)
G_tx = 20;               % Receiver noise (dB)
k = 2*pi/lambda;         % wavenumber
j = sqrt(-1);            % i or j
L1 = 20e-3;                  % length of the center of the rotor
L2 = 100e-3;                  % length of the blade + L1
rpm = 2e4;               % propeller rotations per min
rps = rpm/60*2*pi;       % rotation in rad per sec.

max_vel = rpm/60*2*pi*L2;        % max velocity at rpm
max_dop_shift = 2*max_vel*fc/c;

N_scat = 1;              % N_scat scatters on each side of the origin for each
                         % propeller
N_prop = 4;              % number of propellers
N_blades = 4;            % number of blades
tx_pos = [0; 0; 0];        % Transmitter at origin
tx_vel = [0; 0; 0];        % Transmitter stationary
RCS_bsp = 0.05;          % RCS of blade scattering points
RCS_static = 0.5;        % RCS of static portion of the drone

% contruct a uniformely distributed variable from L1 to L2 or only L2 if
% N_scat = 1
if N_scat > 1
    % bsp := blade scattering position
    bsp = L1:(L2-L1)/(N_scat-1):L2;
else
    bsp = L2;
end

% RCS init; adding one to the end for the static RCS
rcs = zeros(1, N_blades*N_scat*N_prop+1);
rcs(1, :) = RCS_bsp;
rcs(1, length(rcs)) = RCS_static;
% antenna is isotropic so angle doesn't matter. this is only used
% to indicate the number of reflectors.
scat_ang = zeros(1, N_blades*N_scat*N_prop+1);

%% Create Transmitter
tx = phased.Transmitter('PeakPower', txPower, 'Gain', G_tx);

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
rx = phased.ReceiverPreamp('Gain', G_rx, 'NoiseFigure', noiseFigure);

%% Generate CW Signal

cw_waveform_matrix = zeros(N, N_s); % column vector for each subcarrier
offset = floor(N_s/2)+1;

for k=-floor(N_s/2):ceil(N_s/2)-1
    cw_waveform_matrix(1:length(t_symbol),k+offset) = exp(j*2*pi*(fc+k/T_symbol)*t_symbol).';
end

cw_waveform = sum(cw_waveform_matrix,2);
rx_signal = zeros(N);
y_exp = zeros(N, num_iter);
y_theo= zeros(N, num_iter);

% theoretical framework for pos and vel
% standardize notation: first underscore is subscript, second is superscript
x_u_o = [R;0;0]; 
x_t_o = [0;0;0];
v_t_o = [0;0;0];
R_u_o = eul2rotm([0 0 0]);
R_pi_u = eul2rotm([-pi/2 -pi/2 0]);
% let's say the 0th and 2nd go counterclockwise and the other two, clockwise
prop_init_offset_angle = [pi/4;3*pi/4;5*pi/4;7*pi/4];
blade_init_offset_angle = [pi/4;3*pi/4;5*pi/4;7*pi/4];
y = 0.5; % vertical distance between the propellers and the center of the drone
x_bsp_o = zeros(3, N_prop*N_blades*N_scat+1);
v_bsp_o = zeros(3, N_prop*N_blades*N_scat+1);
mag_x_pi_u = 1;
vel_fac = 1;
speed = 2000; % m/s
v_u_o = [0;0;speed];
% set static portion of the drone; the last element in the matrix
x_bsp_o(:, N_prop*N_blades*N_scat+1) = [R;0;0];
v_bsp_o(:, N_prop*N_blades*N_scat+1) = v_u_o;
theo_start = round(2*R/c*fs);

for i = 1:num_iter
    % update scatterers' pos and vel
    time = (i-1)/PRF;
    % blade_ang is the angle of the first blade on the p-th propeller
    blade_ang = blade_init_offset_angle + [time*rps;-time*rps;time*rps;-time*rps];
    x_u_o = [R;0;speed*time];
    v_u_o = [0;0;speed];
    idx = 1;
    
    % N_scat = 1
    for p = 1:N_prop
        theta = blade_ang(p);
        x_pi_u = [mag_x_pi_u*cos(prop_init_offset_angle(p));y;...
            mag_x_pi_u*sin(prop_init_offset_angle(p))];

        for b = 0:N_blades-1
            x_s_pi = bsp*[cos(theta);sin(theta);0];
            x_s_u  = R_pi_u*x_s_pi+x_pi_u;
            x_s_o  = R_u_o*x_s_u+x_u_o;

            v_s_pi = vel_fac*bsp*rps*[-sin(theta);cos(theta);0];
            v_s_o = R_u_o*R_pi_u*v_s_pi + v_u_o;

            x_bsp_o(:,idx) = x_s_o;
            v_bsp_o(:,idx) = v_s_o;
            idx = idx + 1;
            theta = theta + 2*pi/N_blades;
        end
        vel_fac = -1*vel_fac;
    end
    x_bsp_o(:, N_prop*N_blades*N_scat+1) = x_u_o;
    v_bsp_o(:, N_prop*N_blades*N_scat+1) = v_u_o;

    % Transmit Signal
    tx_signal = cw_waveform;
    
    % we need the number of scatters time this signal
    tx_signal = txant(tx_signal, scat_ang);
    prop_signal = channel(tx_signal, x_t_o, x_bsp_o, v_t_o, v_bsp_o);
    refl_signal = target(prop_signal);
    rx_signal= rxant(refl_signal, scat_ang);
    rx_signal = rx(rx_signal);
    y_exp(:,i) = rx_signal;

    % theoretical
    % this is for a theoretical model which uses dedicated scatterers
    idx = 1;
    for t = 0:1/fs:pulse_period
    %idx = theo_start;
    %for t = (2*R/c:1/fs:T_symbol+2*R/c)
        res = 0;
        for k=-floor(N_s/2):ceil(N_s/2)-1
            for p=1:N_prop
                omega = blade_ang(p);
                x_pi_u = [mag_x_pi_u*cos(prop_init_offset_angle(p));y;...
                            mag_x_pi_u*sin(prop_init_offset_angle(p))];
                x_pi_o = R_u_o*x_pi_u + x_u_o;

                for b=1:N_blades
                    x_l2_pi = L2*[cos(omega);sin(omega);0];
                    x_l2_o = R_u_o*R_pi_u*x_l2_pi+x_pi_o;

                    vr = dot(vel_fac*L2*rps*R_u_o*R_pi_u*[-sin(omega);cos(omega);0], x_l2_o/norm(x_l2_o));
                    res = res + RRE(txPower, G_tx, norm(x_pi_o), fc+k/T_symbol, G_rx-noiseFigure)* ... 
                        RCS_bsp*exp(1j*2*pi*(fc+k/T_symbol)*(1+2*vr/c)*(t-2*norm(x_l2_o)/c));
                    omega = omega + 2*pi/N_blades;
                end
                vel_fac = vel_fac*-1;
            end
        end
        res = res + RCS_static*RRE(txPower, G_tx, norm(x_pi_o), fc+k/T_symbol, G_rx-noiseFigure)* ...
            exp(1j*2*pi*(fc+k/T_symbol)*(t-2*norm(x_pi_o)/c));
        y_theo(idx, i) = res;
        idx = idx + 1;
    end
    %}
end



figure('Renderer', 'painters', 'Position', [10 10 900 600])
rdresp  = phased.RangeDopplerResponse('PropagationSpeed',c,'SampleRate',fs,...
    'DopplerFFTLengthSource','Property','DopplerFFTLength',128,'DopplerOutput','Speed',...
    'OperatingFrequency',fc);
mf = cw_waveform(length(t_symbol):-1:1);
mfcoeff = conj(mf);

plotResponse(rdresp,y_exp,mfcoeff);

figure('Renderer', 'painters', 'Position', [10 10 900 600])
mf  = phased.MatchedFilter('Coefficients',mfcoeff);
ymf = mf(y_exp);
[~,ridx] = max(sum(abs(ymf),2)); % detection via peak finding along range
pspectrum(ymf(ridx,:),PRF,'spectrogram')

figure('Renderer', 'painters', 'Position', [10 10 900 600])
plotResponse(rdresp,y_theo,mfcoeff);

figure('Renderer', 'painters', 'Position', [10 10 900 600])
mf  = phased.MatchedFilter('Coefficients',mfcoeff);
ymf = mf(y_theo);
[~,ridx] = max(sum(abs(ymf),2)); % detection via peak finding along range
pspectrum(y_theo(ridx,:),PRF,'spectrogram')

%{
mf  = phased.MatchedFilter('Coefficients',mfcoeff(1:128));
ymf = mf(y_exp);
[~,ridx] = max(sum(abs(ymf),2)); % detection via peak finding along range
pspectrum(ymf(ridx,:),PRF,'spectrogram')
%}

%wav = phased.RectangularWaveform('SampleRate',fs,'PulseWidth',T_symbol,'PRF',PRF);
%{
rdresp  = phased.RangeDopplerResponse('PropagationSpeed',c,'SampleRate',fs,...
    'DopplerFFTLengthSource','Property','DopplerFFTLength',128,'DopplerOutput','Speed',...
    'OperatingFrequency',fc);
mfcoeff = zeros(20, 1);
plotResponse(rdresp,y(:,1:128),mfcoeff);
ylim([0 3000])
%}
% generate matched filter in 
%mfcoeff = conj(cw_waveform(length(cw_waveform):1));

%plot(1:length(cw_waveform),y_exp(:, 1));
%plot(1:length(cw_waveform),cw_waveform);
%x=1;

%{

% theoretical:
theo_sig_matrix = zeros(length(t), N_s);


idx = 1;
vel_fac = 1;
for time=t
    blade_ang = blade_init_offset_angle + [time*rps;-time*rps;time*rps;-time*rps];
    res = 0;

    for k=-floor(N_s/2):ceil(N_s/2)-1
        for p=1:N_prop
            x_pi_u = [mag_x_pi_u*cos(prop_init_offset_angle(p));y;...
                        mag_x_pi_u*sin(prop_init_offset_angle(p))];
            x_pi_o = R_u_o*x_pi_u + x_u_o;

            % sum over n 
            sn = 0;
            for b=1:N_blades
                omega = blade_ang(b);
                q = vel_fac*rps*R_u_o*R_pi_u*[-sin(omega);cos(omega);0];
                l = 4*pi*(fc+k/T_symbol)/c*(dot(q, x_pi_o/norm(x_pi_o)))*(time-2*norm(x_pi_o)/c);
                
                sn = sn + exp(1j*(L1+L2)*l/2)*sinc((L2-L1)*l/2);
            end
    
            res = res + sn*exp(-1j*4*pi*(fc+k/T_symbol)*norm(x_pi_o)/c) + ...
                theo_rcs_ratio*RCS_bsp/RCS_static*exp(1j*2*pi*(fc+k/T_symbol)*(time-2*R/c));
            vel_fac = vel_fac*-1;
        end
    end

    theo_sig_matrix(idx, k+offset) = res;
    idx = idx + 1;
end

theo_sig = sum(theo_sig_matrix, 2);

% shift to baseband
rx_signal_bb = rx_signal .* exp(-1j*2*pi*fc*t).';
cw_waveform_bb = cw_waveform .* exp(-1j*2*pi*fc*t).';
theo_sig_bb = theo_sig .* exp(-1j*2*pi*fc*t).';

%plot time domain

figure;
subplot(3,1,1);
plot(t, cw_waveform_bb, 'b');

subplot(3,1,2);
plot(t, rx_signal_bb, 'r');

subplot(3,1,3);
plot(t, theo_sig_bb, 'g');


% Compute and plot spectrum
%{
start_idx = ceil(2*(R+L2)/c * 1e6 * length(t));
end_idx = length(t);
N_samples = end_idx - start_idx;
%}

% relevant section; only take the portion of the signal from the reflection
t_relevant = t_symbol+2*R/c;

N = length(t_relevant);
T = 1/fs; % time (s)/sample
n_relevant = floor(t_relevant/T);

TX_FFT = abs(fftshift(fft(cw_waveform_bb(n_relevant), N)));
RX_FFT = abs(fftshift(fft(rx_signal_bb(n_relevant), N)));

TX_FFT = TX_FFT/max(TX_FFT);
RX_FFT = RX_FFT/max(RX_FFT);

% theoretical
rx_fft_theo = abs(fftshift(fft(theo_sig_bb(n_relevant), N)));
rx_fft_theo = rx_fft_theo/max(rx_fft_theo);

fft_spacing = -fs/2+1:fs/N:fs/2;

figure('Renderer', 'painters', 'Position', [10 10 900 600])
%plot(f/1e6, TX_FFT, 'b');
plot(fft_spacing/1e6, rx_fft_theo, 'g', fft_spacing/1e6, RX_FFT, 'r', fft_spacing/1e6, TX_FFT, 'b');
xline(-2*max_vel*fc/c/1e6);
xline(2*max_vel*fc/c/1e6);
legend('Rx theoretical','Rx experimental', 'Transmitted',  'Max + velocity', 'Max - velocity');
xlabel('Frequency (MHz)'); ylabel('Normalized Magnitude');
title('Baseband Frequency Spectrum');

%{
% second plot where each is separated
figure('Renderer', 'painters', 'Position', [10 10 900 600])
subplot(3,1,1);
plot(fft_spacing/1e6, rx_fft_theo, 'g');

subplot(3,1,2);
plot(fft_spacing/1e6, RX_FFT, 'r');

subplot(3,1,3);
plot(fft_spacing/1e6, TX_FFT, 'b');
%}



% show STFT with window = 7 us
%{
figure;
stft(rx_signal,fs, FFTLength=2^12);
%spectrogram(rx_signal);
%}

%xlabel('Frequency (MHz)'); ylabel('Magnitude');
%title('Theoretical Frequency Spectrum (Doppler Shift Observed)');
%xlim([900 1100]);
%}


function rx_ampl = RRE(Tx_power, Tx_gain, range, fc, Rx_gain)
    tx_gain_lin = 10^(Tx_gain/10);
    fspl = (4*pi*(2*range)*fc/3e8)^2;
    rx_gain_lin = 10^(Rx_gain/10);
    rx_ampl = Tx_power*tx_gain_lin/fspl*rx_gain_lin;
end