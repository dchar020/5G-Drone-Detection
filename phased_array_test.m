clear all;

% Tutorial:
% https://www.mathworks.com/help/phased/gs/end-to-end-radar-system.html

% generate the waveform
waveform = phased.RectangularWaveform('PulseWidth',1e-6, ...
    'PRF',5e3,'OutputFormat','Pulses','NumPulses',1);

% display the waveform
%Y = waveform();
%plot(Y)

% define the antenna to use
antenna = phased.IsotropicAntennaElement('FrequencyRange',[1e9 10e9]);

% define the target, propagation characteristics and operating freq
target = phased.RadarTarget('Model','Nonfluctuating','MeanRCS',0.5, ...
    'PropagationSpeed',physconst('LightSpeed'),'OperatingFrequency',4e9);

% define the position and velocities of the antenna and target
antennaplatform = phased.Platform('InitialPosition',[0;0;0],'Velocity',[0;0;0]);
targetplatform = phased.Platform('InitialPosition',[7000; 5000; 0], ...
    'Velocity',[-15;-10;0]);

% get the target range and target angles (azim and elev)
[tgtrng,tgtang] = rangeangle(targetplatform.InitialPosition, ...
    antennaplatform.InitialPosition);

% modelling the transmission
Pd = 0.9;
Pfa = 1e-6;
numpulses = 10;
SNR = albersheim(Pd,Pfa,numpulses);

% SNR here is the required SNR for those params (probability of detection
% and probability of false alarm

% radar range equation for determining the required tx power
maxrange = 1.5e4;
lambda = physconst('LightSpeed')/target.OperatingFrequency;
tau = waveform.PulseWidth;
Ts = 290;
rcs = 0.5;
Gain = 20;
dbterm = db2pow(SNR - 2*Gain);
Pt = (4*pi)^3*physconst('Boltzmann')*Ts/tau/rcs/lambda^2*maxrange^4*dbterm;

% Pt is roughly 45 kW. use 50 kW to be safe. Model the tranmission
% using the phase array object
transmitter = phased.Transmitter('PeakPower',50e3,'Gain',20,'LossFactor',0, ...
    'InUseOutputPort',true,'CoherentOnTransmit',true);

% model the radiator and collector;
% radiator: "A radiator converts signals into radiated wavefields"
% collector is the opposite
radiator = phased.Radiator('Sensor',antenna,...
    'PropagationSpeed',physconst('LightSpeed'),'OperatingFrequency',4e9);
collector = phased.Collector('Sensor',antenna,...
    'PropagationSpeed',physconst('LightSpeed'),'Wavefront','Plane', ...
    'OperatingFrequency',4e9);

% model the receiver:
receiver = phased.ReceiverPreamp('Gain',20,'NoiseFigure',2, ...
    'ReferenceTemperature',290,'SampleRate',1e6, ...
    'EnableInputPort',true,'SeedSource','Property','Seed',1e3);

% model the propagation channel:
channel = phased.FreeSpace(...
    'PropagationSpeed',physconst('LightSpeed'), ...
    'OperatingFrequency',4e9,'TwoWayPropagation',false, ...
    'SampleRate',1e6);

% all the parameters and handles are setup. the following code actions
% them

T = 1/waveform.PRF;
% Get antenna position
txpos = antennaplatform.InitialPosition;
% Allocate array for received echoes
rxsig = zeros(waveform.SampleRate*T,numpulses);

for n = 1:numpulses
    % Update the target position
    [tgtpos,tgtvel] = targetplatform(T);
    % Get the range and angle to the target
    [tgtrng,tgtang] = rangeangle(tgtpos,txpos);
    % Generate the pulse
    sig = waveform();
    % Transmit the pulse. Output transmitter status
    [sig,txstatus] = transmitter(sig);
    % Radiate the pulse toward the target
    sig = radiator(sig,tgtang);
    % Propagate the pulse to the target in free space
    sig = channel(sig,txpos,tgtpos,[0;0;0],tgtvel);
    % Reflect the pulse off the target
    sig = target(sig);
    % Propagate the echo to the antenna in free space
    sig = channel(sig,tgtpos,txpos,tgtvel,[0;0;0]);
    % Collect the echo from the incident angle at the antenna
    sig = collector(sig,tgtang);
    % Receive the echo at the antenna when not transmitting
    rxsig(:,n) = receiver(sig,~txstatus);
end

%{
res = zeros(200, 1);
for i = 1:numpulses
    %subplot(numpulses, 1, i)
    %plot(1:200, abs(rxsig(:, i)))
    res = res + abs(rxsig(:, i));
end
plot (res)
%}

% pulsint is just adding all the resulting signals for each pulse.
% this instruction reformats rxsig from [200, 10] to [200, 1]
rxsig = pulsint(rxsig,'noncoherent');

% unigrid is uniform grid; a generator for basic expressions
t = unigrid(0,1/receiver.SampleRate,T,'[)');

% range gates are the range bins; there's the same amount of bins as there
% are time intervals for this sampling rate
rangegates = (physconst('LightSpeed')*t)/2;

% plot.
plot(rangegates/1e3,rxsig)
hold on
xlabel('range (km)')
ylabel('Power')
xline([tgtrng/1e3,tgtrng/1e3],'r')
hold off
