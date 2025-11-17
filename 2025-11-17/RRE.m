% Parameters
c = physconst('LightSpeed');
fc = 50e9;                     % 10 GHz
lambda = c / fc;
range = 100;                 % meters
rcs = 1;                      % m^2
txPower = 3.14;                % 1000 W = 60 dBm
G_tx = 10;                    % dB
G_rx = 10;                    % dB
L = 1;                        % System loss factor (unitless)
fs = 1e6;                     % 1 MHz sample rate
noiseFigure = 0;              % 0 dB for clean comparison

% Radar Equation (Linear form)
Pr = (txPower * 10^(G_tx/10) * 10^(G_rx/10) * lambda^2 * rcs) ...
     / ((4*pi)^3 * range^4 * L);
Pr_dB = 10*log10(Pr);

% Print theoretical received power
fprintf('Theoretical received power: %.2f dBW\n', Pr_dB);

%% Create System Components
tx = phased.Transmitter('PeakPower', txPower, 'Gain', G_tx);
ant = phased.IsotropicAntennaElement('BackBaffled', true);
radiator = phased.Radiator('Sensor', ant, 'PropagationSpeed', c, 'OperatingFrequency', fc);
collector = phased.Collector('Sensor', ant, 'PropagationSpeed', c, 'OperatingFrequency', fc);
target = phased.RadarTarget('MeanRCS', rcs, 'OperatingFrequency', fc, 'Model', 'Nonfluctuating');
channel = phased.FreeSpace('SampleRate', fs, 'TwoWayPropagation', true, 'OperatingFrequency', fc);
rx = phased.ReceiverPreamp('Gain', G_rx, 'NoiseFigure', noiseFigure);

% Transmit a simple pulse
N = 100;                                 % Pulse length
txPulse = ones(N,1);                     % Flat pulse
txOut = tx(txPulse);
txRadiated = radiator(txOut, [0;0]);

% Propagation
target_pos = [range; 0; 0];              % Target straight ahead
tx_pos = [0; 0; 0];
v_t = [0; 0; 0];
v_r = [0; 0; 0];
prop_sig = channel(txRadiated, tx_pos, target_pos, v_t, v_r);

% Target reflection
tgt_reflected = target(prop_sig);

% Return propagation
%rx_sig = channel(tgt_reflected, target_pos, tx_pos, v_t, v_r);

% Collect and receive
rx_collected = collector(tgt_reflected, [0;0]);
rx_out = rx(rx_collected);

% Estimate power (in dBW)
rx_meas_dB = 10*log10(mean(abs(rx_out).^2));
fprintf('Simulated received power:  %.2f dBW\n', rx_meas_dB);

% Difference
fprintf('Difference: %.2f dB\n', rx_meas_dB - Pr_dB);
