clc; clear; close all;

j = sqrt(-1);            % i or j

% Electromagnetic parameters
fc = 50e9;               % Carrier frequency (50GHz)
c = 3e8;                 % Speed of light (m/s)
lambda = c/fc;           % Wavelength
T_symbol = 1/(960e3);    % OFDM symbol duration for 5G. Smallest possible value
N_s = 14;                % num subcarriers
SCS = 1/T_symbol;        % subcarrier spacing
BW = SCS*(N_s+1);	     % Bandwidth
k = 2*pi/lambda;         % wavenumber

% Discretization parameters
fs = BW;                 % Sampling frequency (narrowband sampling < Nyquist)
N_sps = T_symbol*fs;	 % Number of Samples Per OFDM Symbol
M = 2;					 % Number of symbols in one pulse
N = M*N_sps;	         % Number of samples in one pulse 


% Physical parameters 
R = 100;                        % Target range (100 m)
v = 0;                          % Target velocity (0 m/s)
L1 = 20e-3;                     % length of the center of the rotor
L2 = 100e-3;                    % length of the blade + L1
rpm = 2e4;                      % propeller rotations per min
rps = rpm/60*2*pi;              % rotation in rad per sec.
max_vel = rpm/60*2*pi*L2;       % max velocity at rpm
max_dop_shift = 2*max_vel*fc/c; % max doppler shift
N_scat = 1;              % N_scat scatters on each side of the origin for each
                         % propeller
N_prop = 4;              % number of propellers
N_blades = 4;            % number of blades
tx_pos = [0; 0; 0];      % Transmitter at origin
tx_vel = [0; 0; 0];      % Transmitter stationary
RCS_bsp = 0.05;          % RCS of blade scattering points
RCS_static = 0.1;        % RCS of static portion of the drone
blade_width = 0.004;     % width of the propeller blade

% Radar parameters 
pulse_period = T_symbol*M;
PRF = 1/pulse_period;
txPower = 3.16;          % Transmitted power (Watts)
noiseFigure = 6;         % Receiver noise figure (dB)
G_rx = 10;               % Receiver gain (dBi)
G_tx = 10;               % Receiver noise (dBi)

% Simulation parameters 
num_iter = 2e3; % number of pulses
t_symbol = 0:1/fs:T_symbol; % interval for the first symbol

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

% Full-system simulation using Phased Array
%-------------------------------------------------------------------------%
% Create Transmitter
tx = phased.Transmitter('PeakPower', txPower, 'Gain', G_tx);

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
rx = phased.ReceiverPreamp('Gain', G_rx, 'NoiseFigure', noiseFigure);

% Generate CW Signal
cw_waveform_matrix = zeros(N, N_s); % column vector for each subcarrier
offset = floor(N_s/2)+1;

for k=-floor(N_s/2):ceil(N_s/2)-1
    cw_waveform_matrix(1:length(t_symbol),k+offset) = exp(j*2*pi*(fc+k/T_symbol)*t_symbol).';
end

cw_waveform = sum(cw_waveform_matrix,2);
%-------------------------------------------------------------------------%

% Init data storage variables
rx_signal = zeros(N);
y_exp = zeros(N, num_iter);
y_theo= zeros(N, num_iter);
y_theo_ps= zeros(N, num_iter);

% theoretical framework for pos and vel
% standardize notation: first underscore is subscript, second is superscript
x_u_o = [R;0;0]; 
x_t_o = [0;0;0];
v_t_o = [0;0;0];
R_u_o = eul2rotm([0 0 0]);
R_pi_u = eul2rotm([-pi/2 -pi/2 0]);
% the 0th and 2nd go counterclockwise and the other two, clockwise
prop_init_offset_angle = [pi/4;3*pi/4;5*pi/4;7*pi/4];
blade_init_offset_angle = [pi/4;3*pi/4;5*pi/4;7*pi/4];
y = 0.5; % vertical distance between the propellers and the center of the drone
x_bsp_o = zeros(3, N_prop*N_blades*N_scat+1);
v_bsp_o = zeros(3, N_prop*N_blades*N_scat+1);
mag_x_pi_u = 1;
vel_fac = 1;
% set static portion of the drone; the last element in the matrix
x_bsp_o(:, N_prop*N_blades*N_scat+1) = [R;0;0];
v_bsp_o(:, N_prop*N_blades*N_scat+1) = [0;0;0];
theo_start = round(2*R/c*fs);

% Main simulation loop
for i = 1:num_iter
	time = (i-1)/PRF;
    % update scatterers' pos and vel
%-------------------------------------------------------------------------%
    % blade_ang is the angle of the first blade on the p-th propeller
    blade_ang = blade_init_offset_angle + [time*rps;-time*rps;time*rps;-time*rps];
    idx = 1;
    
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

                    v_s_pi = vel_fac*bsp*rps*[-sin(theta);cos(theta);0];
                    v_s_o = R_u_o*R_pi_u*v_s_pi;

                    x_bsp_o(:,idx) = x_s_o;
                    v_bsp_o(:,idx) = v_s_o;
                    idx = idx + 1;
                end
            end
            % move to next propeller, reverse direction
            vel_fac = -1*vel_fac;
        end
    else
        for p = 1:N_prop
            theta = blade_ang(p);
            x_pi_u = [mag_x_pi_u*cos(prop_init_offset_angle(p));y;...
                mag_x_pi_u*sin(prop_init_offset_angle(p))];

            for b = 0:N_blades-1
                x_s_pi = bsp*[cos(theta);sin(theta);0];
                x_s_u  = R_pi_u*x_s_pi+x_pi_u;
                x_s_o  = R_u_o*x_s_u+x_u_o;

                v_s_pi = vel_fac*bsp*rps*[-sin(theta);cos(theta);0];
                v_s_o = R_u_o*R_pi_u*v_s_pi;

                x_bsp_o(:,idx) = x_s_o;
                v_bsp_o(:,idx) = v_s_o;
                idx = idx + 1;
                theta = theta + 2*pi/N_blades;
            end
            vel_fac = -1*vel_fac;
        end
    end
%-------------------------------------------------------------------------%
	% Phased Array package processing:
    % Transmit Signal
    tx_signal = cw_waveform;
    
    % we need the number of scatters time this signal
    tx_signal = txant(tx(tx_signal), scat_ang);
    prop_signal = channel(tx_signal, x_t_o, x_bsp_o, v_t_o, v_bsp_o);
    refl_signal = target(prop_signal, true);
    rx_signal= rxant(refl_signal, scat_ang);
    rx_signal = rx(rx_signal);
    y_exp(:,i) = rx_signal;

    % theoretical MM model
    idx = 1;
    for t = (0:1/fs:pulse_period)
        res = 0;
		
        for k=-floor(N_s/2):ceil(N_s/2)-1
            for p=1:N_prop
                omega = -blade_ang(p);
                x_pi_u = [mag_x_pi_u*cos(prop_init_offset_angle(p));y;...
                            mag_x_pi_u*sin(prop_init_offset_angle(p))];
                x_pi_o = R_u_o*x_pi_u + x_u_o;
    
                % sum over n 
                sn = 0;
                for b=1:N_blades
                    x_l2_pi = L2*[cos(omega);sin(omega);0];
                    x_l2_o = R_u_o*R_pi_u*x_l2_pi+x_pi_o;

                    q = vel_fac*rps*R_u_o*R_pi_u*[-sin(omega);cos(omega);0];
                    l = 4*pi*(fc+k/T_symbol)/c*(dot(q, x_l2_o/norm(x_l2_o)))*(t-2*norm(x_l2_o)/c);
                    
                    res = res +  (L2-L1)*RRE(txPower, G_tx, norm(x_l2_o), fc+k/T_symbol, G_rx, 1+blade_width)*...
                        exp(1j*2*pi*(fc+k/T_symbol)*(t-2*norm(x_l2_o)/c))*exp(1j*(L1+L2)*l/2)*sinc((L2-L1)*l/2);
                    omega = omega + 2*pi/N_blades;
                end
                res = res + RRE(txPower, G_tx, norm(x_u_o), fc+k/T_symbol, G_rx, RCS_static)*...
                    exp(1j*2*pi*(fc+k/T_symbol)*(t-2*R/c));
                vel_fac = vel_fac*-1;
            end
        end
        
        if (2*(R-5)/c <= t) & (t <= T_symbol+2*(R+5)/c)
            y_theo(idx, i) = res/N_s;
        end
        idx = idx + 1;
    end

    
    % theoretical, PS model    
    idx = 1;
    for t = 0:1/fs:pulse_period
        res = 0;
        for k=-floor(N_s/2):ceil(N_s/2)-1
            for p=1:N_prop
                omega = -blade_ang(p);
                x_pi_u = [mag_x_pi_u*cos(prop_init_offset_angle(p));y;...
                            mag_x_pi_u*sin(prop_init_offset_angle(p))];
                x_pi_o = R_u_o*x_pi_u + x_u_o;

                for b=1:N_blades
                    x_l2_pi = L2*[cos(omega);sin(omega);0];
                    x_l2_o = R_u_o*R_pi_u*x_l2_pi+x_pi_o;

                    vr = dot(vel_fac*L2*rps*R_u_o*R_pi_u*[-sin(omega);cos(omega);0], x_l2_o/norm(x_l2_o));
                    res = res + RRE(txPower, G_tx, norm(x_l2_o), fc+k/T_symbol, G_rx, RCS_bsp)* ... 
                        exp(1j*2*pi*(fc+k/T_symbol)*(1+2*vr/c)*(t-2*norm(x_l2_o)/c));
                    omega = omega + 2*pi/N_blades;
                end
                vel_fac = vel_fac*-1;
            end
            res = res + RRE(txPower, G_tx, norm(x_u_o), fc+k/T_symbol, G_rx, RCS_static)* ...
             exp(1j*2*pi*(fc+k/T_symbol)*(t-2*norm(x_pi_o)/c));
        end
        
        if (2*(R-5)/c <= t) & (t <= T_symbol+2*(R+5)/c)
            y_theo_ps(idx, i) = res/N_s;
        end
        idx = idx + 1;
    end
end


% Generate figures
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
pspectrum(ymf(ridx,:),PRF,'spectrogram')

figure('Renderer', 'painters', 'Position', [10 10 900 600])
plotResponse(rdresp,y_theo_ps,mfcoeff);

figure('Renderer', 'painters', 'Position', [10 10 900 600])
mf  = phased.MatchedFilter('Coefficients',mfcoeff);
ymf = mf(y_theo_ps);
[~,ridx] = max(sum(abs(ymf),2)); % detection via peak finding along range
pspectrum(ymf(ridx,:),PRF,'spectrogram')

% save data
%{
r_fss = y_exp;
r_mm = y_theo;
r_ps = y_theo_ps;

writematrix(r_mm,'DASC_Charron_r_mm.csv');
writematrix(r_ps,'DASC_Charron_r_ps.csv');
writematrix(r_fss,'DASC_Charron_r_fss.csv');
save('DASC_Charron_receive_data.mat', 'r_fss', 'r_mm', 'r_ps');
%}

% Function RRE(): returns an amplitude for the backscatter signal 
function rx_ampl = RRE(Tx_power, Tx_gain, range, fc, Rx_gain, rcs)

    lambda = 3e8/fc;
    L = 1;

    rx_power_ampl = (Tx_power*10^(Tx_gain/10)*10^(Rx_gain/10)*lambda^2*rcs ...
        /((4*pi)^3*range^4*L));

    rx_ampl = sqrt(rx_power_ampl*50);

end

% Additional figure generation and debugging tools if one is interested

% figure;
% plot((1:length(cw_waveform)), cw_waveform);
% title('CW 5G OFDM Transmit Signal Pulse')
% xlabel('Discrete-time Sample Number'); ylabel('Amplitude')
% 
% num_steps = N - length(t_symbol);
% mf  = phased.MatchedFilter('Coefficients',mfcoeff);
% ymf = mf(y_theo);
% y = sum(abs(ymf),2);
% r = c/(2*BW);
% ax = (0:num_steps)*r;
% y_ss = (y((length(y)-num_steps):length(y)));
% th = num_steps*r+r/2;
% a = r/2:r:th;
% 
% figure('Renderer', 'painters', 'Position', [10 10 900 600])
% subplot(3,1,1)
% plot(ax, y_ss);
% grid on
% xticks(a);
% title('MM Model Range Bin Correlation with Transmit Signal')
% xlabel('Range (m)'); ylabel('Correlation')
% 
% subplot(3,1,2)
% ymf = mf(y_theo_ps);
% y = sum(abs(ymf),2);
% y_ss = y((length(y)-num_steps):length(y));
% plot(ax, y_ss);
% grid on
% xticks(a);
% title('PS Model Range Bin Correlation with Transmit Signal')
% xlabel('Range (m)'); ylabel('Correlation')
% 
% subplot(3,1,3)
% ymf = mf(y_exp);
% y = sum(abs(ymf),2);
% y_ss = y((length(y)-num_steps):length(y));
% plot(ax, y_ss);
% grid on
% xticks(a);
% title('Exp Model Range Bin Correlation with Transmit Signal')
% xlabel('Range (m)'); ylabel('Correlation')
% 
% 
% %}
% 
% idx = 2;
% 
% %figure('Renderer', 'painters', 'Position', [10 10 900 600])
% figure('Renderer', 'painters', 'Position', [10 10 1000 900])
% 
% subplot(3,3,1);
% plot(1:length(y_theo(:,idx)), y_theo(:,idx-1));
% title('MM model \textbf{r}$_x[:,1]$', Interpreter='latex')
% subplot(3,3,2);
% plot(1:length(y_theo(:,idx)), y_theo(:,idx));
% title('MM model \textbf{r}$_x[:,2]$', Interpreter='latex')
% subplot(3,3,3);
% plot(1:length(y_theo(:,idx)), y_theo(:,idx+1));
% title('MM model \textbf{r}$_x[:,3]$', Interpreter='latex')
% 
% subplot(3,3,4);
% plot(1:length(y_theo(:,idx)), y_theo_ps(:,idx-1));
% title('PS model \textbf{r}$_x[:,1]$', Interpreter='latex')
% subplot(3,3,5);
% plot(1:length(y_theo(:,idx)), y_theo_ps(:,idx));
% title('PS model \textbf{r}$_x[:,2]$', Interpreter='latex')
% subplot(3,3,6);
% plot(1:length(y_theo(:,idx)), y_theo_ps(:,idx+1));
% title('PS model \textbf{r}$_x[:,3]$', Interpreter='latex')
% 
% subplot(3,3,7);
% plot(1:length(y_exp(:,idx)), y_exp(:,idx-1));
% title('Exp model \textbf{r}$_x[:,1]$', Interpreter='latex')
% xlabel('Discrete-time Sample Number');
% subplot(3,3,8);
% plot(1:length(y_exp(:,idx)), y_exp(:,idx));
% title('Exp model \textbf{r}$_x[:,2]$', Interpreter='latex')
% xlabel('Discrete-time Sample Number');
% subplot(3,3,9);
% plot(1:length(y_exp(:,idx)), y_exp(:,idx+1));
% title('Exp model \textbf{r}$_x[:,3]$', Interpreter='latex')
% xlabel('Discrete-time Sample Number');
% 
% %}
% %{
% mf  = phased.MatchedFilter('Coefficients',mfcoeff);
% ymf = mf(y_theo);
% y = sum(abs(ymf),2);
% [~,ridx] = max(sum(abs(ymf),2)); % detection via peak finding along range
% %}
% 
% %{
% figure('Renderer', 'painters', 'Position', [10 10 900 600])
% plot(1:length(y_theo(ridx,:)), y_theo(ridx,:)/abs(max(y_theo(ridx,:))));
% xlabel('Discrete Sample Number'); ylabel('Normalized Magnitude');
% title('Theoretical MM Time-Domain Radar Signal');
% %}
% %{
% figure('Renderer', 'painters', 'Position', [10 10 900 600])
% plotResponse(rdresp,y_theo,mfcoeff);
% %}
% %{
% figure;
% mf  = phased.MatchedFilter('Coefficients',mfcoeff);
% ymf = mf(y_theo);
% [~,ridx] = max(sum(abs(ymf),2)); % detection via peak finding along range
% pspectrum(y_theo(ridx-1,:),PRF,'spectrogram')
% %}
% 
% %{
% M = 80;
% L = 16;
% lk = 0.7;
% x = abs(fftshift(fft(ymf(ridx-1,1:M), 80)));
% 
% figure('Renderer', 'painters', 'Position', [10 10 900 600])
% plot(1:length(x), x);
% %}
% 
% 
% figure('Renderer', 'painters', 'Position', [10 10 900 600])
% subplot(1,3,1)
% mf  = phased.MatchedFilter('Coefficients',mfcoeff);
% ymf = mf(y_theo);
% [~,ridx] = max(sum(abs(ymf),2)); % detection via peak finding along range
% pspectrum(ymf(ridx,:),PRF,'spectrogram');
% 
% subplot(1,3,2)
% mf  = phased.MatchedFilter('Coefficients',mfcoeff);
% ymf = mf(y_theo_ps);
% [~,ridx] = max(sum(abs(ymf),2)); % detection via peak finding along range
% pspectrum(ymf(ridx,:),PRF,'spectrogram');
% 
% subplot(1,3,3)
% mf  = phased.MatchedFilter('Coefficients',mfcoeff);
% ymf = mf(y_exp);
% [~,ridx] = max(sum(abs(ymf),2)); % detection via peak finding along range
% pspectrum(ymf(ridx,:),PRF,'spectrogram');
% 
% rdresp  = phased.RangeDopplerResponse('PropagationSpeed',c,'SampleRate',fs,...
%     'DopplerFFTLengthSource','Property','DopplerFFTLength',128,'DopplerOutput','Speed',...
%     'OperatingFrequency',fc);
% mf = cw_waveform(length(t_symbol):-1:1);
% mfcoeff = conj(mf);
% 
% figure('Renderer', 'painters', 'Position', [10 10 900 600])
% subplot(1,3,1)
% plotResponse(rdresp,y_theo,mfcoeff);
% %colorbar('Ticks',[-150, -30]);
% 
% subplot(1,3,2)
% plotResponse(rdresp,y_theo_ps,mfcoeff);
% 
% subplot(1,3,3)
% plotResponse(rdresp,y_exp,mfcoeff);